pub const MIN_ORF: usize = 90;
pub const N_HEX: usize = 4096; // 4^6
pub const N_TRI: usize = 64;   // 4^3

pub type HexModel = [f64; N_HEX];
pub type MonoModel = [f64; N_TRI];

#[derive(Clone)]
pub struct Gene {
    pub start: usize,        // 1-based genome coord (smaller)
    pub end: usize,          // 1-based genome coord (larger)
    pub is_plus: bool,
    pub seq_start: usize,    // 0-based on strand sequence
    pub seq_end: usize,      // 0-based exclusive
    pub length: usize,
    pub stop_group: u32,
    pub frame: u8,
    pub start_codon: [u8; 3],
    pub is_longest: bool,
    // Scores
    pub hex_avg: f64,
    pub hex_total: f64,
    pub frame_bias: f64,
    pub hex_cov: f64,
    pub rbs: f64,
    pub rbs_pwm: f64,
    pub gc3: f64,
    pub upstream_at: f64,
    pub mono: f64,
    pub edge: f64,
    pub score: f64,
    pub weight: f64,
    pub shadow_pen: f64,
    pub _sq: f64,
    pub _base_weight: f64,
    pub adj_bonus: f64,
    pub leaderless: f64,  // -10 box score (cycle 20 idea)
    pub start_ctx: f64,   // start codon context score (cycle 25 idea)
    pub gc3_bias: f64,    // |GC3 - GC12| wobble position bias (cycle 15 idea)
    pub viterbi_frac: f64, // fraction of region marked coding by nucleotide-level DP (new signal)
    pub start_nn: f64,     // start context neural model score
    pub stop_nn: f64,      // stop context neural model score
    pub upstream_gene_dist: i64, // distance to nearest upstream gene (same strand), -1 if overlap
    pub upstream_coding: f64,    // coding potential of upstream in-frame region (truncation evidence)
}

impl Gene {
    pub fn new() -> Self {
        Gene {
            start: 0, end: 0, is_plus: true,
            seq_start: 0, seq_end: 0, length: 0,
            stop_group: 0, frame: 0,
            start_codon: [0; 3], is_longest: false,
            hex_avg: 0.0, hex_total: 0.0,
            frame_bias: 0.0, hex_cov: 0.5,
            rbs: 0.0, rbs_pwm: 0.0,
            gc3: 0.5, upstream_at: 0.5,
            mono: 0.0,
            edge: 0.0, score: 0.0, weight: 0.0,
            shadow_pen: 1.0,
            _sq: 0.0, _base_weight: 0.0, adj_bonus: 0.0, leaderless: 0.0, start_ctx: 0.0, gc3_bias: 0.0, viterbi_frac: 0.0, start_nn: 0.0, stop_nn: 0.0, upstream_gene_dist: i64::MAX, upstream_coding: 0.0,
        }
    }

    pub fn start_type(&self) -> f64 {
        match &self.start_codon {
            b"ATG" => 1.0,
            b"GTG" => 0.55,
            b"TTG" => 0.35,
            _ => 0.2,
        }
    }

    pub fn is_atg(&self) -> bool { self.start_codon == *b"ATG" }
}
