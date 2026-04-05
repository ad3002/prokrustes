//! LightGBM model inference in pure Rust.
//! Parses the text model file and evaluates decision trees.
//! No external dependencies.

use std::fs;

/// A single decision tree node.
struct TreeNode {
    split_feature: usize,
    threshold: f64,
    left_child: i32,  // negative = leaf index (-(idx+1)), positive = node index
    right_child: i32,
    leaf_value: f64,   // only used if this is a leaf
}

/// A single decision tree.
struct Tree {
    nodes: Vec<TreeNode>,
    leaves: Vec<f64>,
    num_leaves: usize,
}

impl Tree {
    fn predict(&self, features: &[f64]) -> f64 {
        let mut node_idx: i32 = 0;
        loop {
            if node_idx < 0 {
                // Leaf node
                let leaf_idx = (-node_idx - 1) as usize;
                return if leaf_idx < self.leaves.len() { self.leaves[leaf_idx] } else { 0.0 };
            }
            let ni = node_idx as usize;
            if ni >= self.nodes.len() { return 0.0; }
            let node = &self.nodes[ni];
            let feat_val = if node.split_feature < features.len() {
                features[node.split_feature]
            } else {
                0.0
            };
            // decision_type=2 means <= (LightGBM default for numerical)
            node_idx = if feat_val <= node.threshold {
                node.left_child
            } else {
                node.right_child
            };
        }
    }
}

/// LightGBM ensemble model.
pub struct LgbModel {
    trees: Vec<Tree>,
    pub feature_names: Vec<String>,
    shrinkage: f64,
}

impl LgbModel {
    /// Load model from LightGBM text file.
    pub fn load(path: &str) -> Option<Self> {
        let content = fs::read_to_string(path).ok()?;
        let mut feature_names = Vec::new();
        let mut trees = Vec::new();
        let mut global_shrinkage = 0.05;

        let lines: Vec<&str> = content.lines().collect();
        let mut i = 0;

        while i < lines.len() {
            let line = lines[i].trim();

            if line.starts_with("feature_names=") {
                feature_names = line[14..].split_whitespace()
                    .map(|s| s.to_string()).collect();
            }

            if line.starts_with("Tree=") {
                // Parse one tree
                let mut num_leaves = 0;
                let mut split_feature = Vec::new();
                let mut threshold = Vec::new();
                let mut left_child = Vec::new();
                let mut right_child = Vec::new();
                let mut leaf_value = Vec::new();
                let mut shrinkage_val = 0.05;

                i += 1;
                while i < lines.len() {
                    let tl = lines[i].trim();
                    if tl.is_empty() || tl.starts_with("Tree=") || tl == "end of trees" {
                        break;
                    }
                    if tl.starts_with("num_leaves=") {
                        num_leaves = tl[11..].parse().unwrap_or(0);
                    } else if tl.starts_with("split_feature=") {
                        split_feature = tl[14..].split_whitespace()
                            .map(|s| s.parse::<usize>().unwrap_or(0)).collect();
                    } else if tl.starts_with("threshold=") {
                        threshold = tl[10..].split_whitespace()
                            .map(|s| s.parse::<f64>().unwrap_or(0.0)).collect();
                    } else if tl.starts_with("left_child=") {
                        left_child = tl[11..].split_whitespace()
                            .map(|s| s.parse::<i32>().unwrap_or(0)).collect();
                    } else if tl.starts_with("right_child=") {
                        right_child = tl[12..].split_whitespace()
                            .map(|s| s.parse::<i32>().unwrap_or(0)).collect();
                    } else if tl.starts_with("leaf_value=") {
                        leaf_value = tl[11..].split_whitespace()
                            .map(|s| s.parse::<f64>().unwrap_or(0.0)).collect();
                    } else if tl.starts_with("shrinkage=") {
                        shrinkage_val = tl[10..].parse().unwrap_or(0.05);
                        global_shrinkage = shrinkage_val;
                    }
                    i += 1;
                }

                if num_leaves > 0 && !split_feature.is_empty() {
                    let n_internal = split_feature.len();
                    let mut nodes = Vec::with_capacity(n_internal);
                    for j in 0..n_internal {
                        nodes.push(TreeNode {
                            split_feature: split_feature[j],
                            threshold: threshold.get(j).copied().unwrap_or(0.0),
                            left_child: left_child.get(j).copied().unwrap_or(0),
                            right_child: right_child.get(j).copied().unwrap_or(0),
                            leaf_value: 0.0,
                        });
                    }
                    trees.push(Tree { nodes, leaves: leaf_value, num_leaves });
                }
                continue;
            }
            i += 1;
        }

        if trees.is_empty() { return None; }
        eprintln!("LGB model: {} trees, {} features, shrinkage={}",
                  trees.len(), feature_names.len(), global_shrinkage);
        Some(LgbModel { trees, feature_names, shrinkage: global_shrinkage })
    }

    /// Predict raw score for a feature vector.
    pub fn predict(&self, features: &[f64]) -> f64 {
        self.trees.iter().map(|t| t.predict(features)).sum()
    }
}
