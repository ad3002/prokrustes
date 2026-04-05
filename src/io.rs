use std::fs;

pub fn read_fasta(path: &str) -> Vec<u8> {
    let data = fs::read_to_string(path).expect("Cannot read FASTA file");
    let mut seq = Vec::with_capacity(5_000_000);
    for line in data.lines() {
        if !line.starts_with('>') {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    seq
}

pub fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N'
    }).collect()
}

#[inline(always)]
pub fn nt4(b: u8) -> usize {
    match b { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4 }
}

#[inline]
pub fn hex_enc(s: &[u8]) -> Option<usize> {
    let mut v = 0usize;
    for i in 0..6 {
        let n = nt4(s[i]);
        if n > 3 { return None; }
        v = v * 4 + n;
    }
    Some(v)
}

#[inline]
pub fn tri_enc(s: &[u8]) -> Option<usize> {
    let mut v = 0usize;
    for i in 0..3 {
        let n = nt4(s[i]);
        if n > 3 { return None; }
        v = v * 4 + n;
    }
    Some(v)
}
