#!/usr/bin/env python3
"""Train LightGBM pairwise start ranker from exported training data.
Uses lambdarank objective: learns to rank starts within stop groups."""

import os
import numpy as np

def main():
    data_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'start_training_v2.tsv')
    model_path = os.path.join(os.path.dirname(__file__), '..', 'models', 'start_ranker_v2.lgb')

    if not os.path.exists(data_path):
        print(f"Training data not found: {data_path}")
        print("Run: python3 scripts/export_start_training.py")
        return

    # Load data
    qids = []
    labels = []
    features = []
    feature_names = None

    with open(data_path) as f:
        header = f.readline().strip().split('\t')
        feature_names = header[2:]  # skip qid, label
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            qids.append(int(parts[0]))
            labels.append(int(parts[1]))
            features.append([float(x) for x in parts[2:]])

    X = np.array(features)
    y = np.array(labels)
    qids_arr = np.array(qids)

    # Group sizes for lambdarank
    unique_qids = sorted(set(qids_arr))
    group_sizes = [int(np.sum(qids_arr == q)) for q in unique_qids]

    print(f"Data: {len(y)} samples, {len(unique_qids)} groups, {X.shape[1]} features")
    print(f"Positive labels: {y.sum()} ({y.mean()*100:.1f}%)")
    print(f"Group sizes: min={min(group_sizes)}, median={sorted(group_sizes)[len(group_sizes)//2]}, max={max(group_sizes)}")

    # Split: hold out 20% of groups for validation
    n_train = int(len(unique_qids) * 0.8)
    train_qids = set(unique_qids[:n_train])
    val_qids = set(unique_qids[n_train:])

    train_mask = np.array([q in train_qids for q in qids_arr])
    val_mask = ~train_mask

    X_train, y_train = X[train_mask], y[train_mask]
    X_val, y_val = X[val_mask], y[val_mask]

    train_groups = [int(np.sum(qids_arr[train_mask] == q)) for q in sorted(train_qids)]
    val_groups = [int(np.sum(qids_arr[val_mask] == q)) for q in sorted(val_qids)]

    print(f"Train: {len(y_train)} samples, {len(train_groups)} groups")
    print(f"Val:   {len(y_val)} samples, {len(val_groups)} groups")

    try:
        import lightgbm as lgb
    except ImportError:
        print("\nLightGBM not installed. Install: pip install lightgbm")
        print("Using simple logistic regression instead...")
        train_simple(X, y, qids_arr, feature_names, model_path)
        return

    # Train LightGBM lambdarank
    train_data = lgb.Dataset(X_train, label=y_train, group=train_groups,
                             feature_name=feature_names)
    val_data = lgb.Dataset(X_val, label=y_val, group=val_groups,
                           feature_name=feature_names, reference=train_data)

    params = {
        'objective': 'lambdarank',
        'metric': 'ndcg',
        'ndcg_eval_at': [1],
        'num_leaves': 31,
        'learning_rate': 0.05,
        'min_data_in_leaf': 5,
        'feature_fraction': 0.8,
        'verbose': -1,
    }

    model = lgb.train(params, train_data, num_boost_round=300,
                      valid_sets=[val_data],
                      callbacks=[lgb.early_stopping(20), lgb.log_evaluation(50)])

    # Evaluate
    val_pred = model.predict(X_val)
    # Per-group accuracy: does the highest-scored candidate have label=1?
    val_qids_sorted = sorted(val_qids)
    correct = 0
    total = 0
    for q in val_qids_sorted:
        mask = qids_arr == q
        if not mask.any(): continue
        q_labels = y[mask]
        q_preds = model.predict(X[mask])
        if q_labels.sum() == 0: continue  # no correct start in this group
        best_pred_idx = np.argmax(q_preds)
        if q_labels[best_pred_idx] == 1:
            correct += 1
        total += 1

    print(f"\nValidation: top-1 accuracy = {correct}/{total} ({correct*100//max(total,1)}%)")

    # Feature importance
    importance = model.feature_importance(importance_type='gain')
    sorted_idx = np.argsort(-importance)
    print(f"\nTop features:")
    for idx in sorted_idx[:10]:
        print(f"  {feature_names[idx]:>25}: {importance[idx]:.0f}")

    # Save model
    model.save_model(model_path, num_iteration=model.best_iteration)
    print(f"\nModel saved: {model_path}")

    # Also save as text for Rust parser
    text_path = model_path.replace('.lgb', '.txt')
    model.save_model(text_path, num_iteration=model.best_iteration)
    print(f"Text model: {text_path}")

def train_simple(X, y, qids, feature_names, model_path):
    """Fallback: simple feature analysis without LightGBM."""
    import statistics

    print("\n=== Feature importance (simple analysis) ===")
    unique_qids = sorted(set(qids))

    for fi, fname in enumerate(feature_names):
        # For each group: does correct start have higher value?
        correct_higher = 0
        total = 0
        for q in unique_qids:
            mask = qids == q
            q_labels = y[mask]
            q_feat = X[mask, fi]
            if q_labels.sum() != 1: continue
            correct_idx = q_labels.argmax()
            correct_val = q_feat[correct_idx]
            best_other = max(q_feat[j] for j in range(len(q_feat)) if j != correct_idx) if len(q_feat) > 1 else correct_val
            if correct_val > best_other + 0.001:
                correct_higher += 1
            total += 1

        if total > 0:
            pct = correct_higher * 100 // total
            print(f"  {fname:>25}: correct > best_other in {correct_higher}/{total} ({pct}%)")

if __name__ == "__main__":
    main()
