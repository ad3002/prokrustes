use std::fs;

/// Read multi-FASTA file. Returns vector of (name, sequence) pairs.
pub fn read_fasta_multi(path: &str) -> Vec<(String, Vec<u8>)> {
    let data = fs::read_to_string(path).expect("Cannot read FASTA file");
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::with_capacity(5_000_000);

    for line in data.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !name.is_empty() || !seq.is_empty() {
                records.push((name.clone(), std::mem::take(&mut seq)));
            }
            name = header.split_whitespace().next().unwrap_or("").to_string();
            seq = Vec::with_capacity(5_000_000);
        } else {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !seq.is_empty() {
        records.push((name, seq));
    }
    records
}

/// Read FASTA file as single concatenated sequence (legacy, for single-contig genomes).
pub fn read_fasta(path: &str) -> Vec<u8> {
    let records = read_fasta_multi(path);
    let total: usize = records.iter().map(|(_, s)| s.len()).sum();
    let mut seq = Vec::with_capacity(total);
    for (_, s) in records {
        seq.extend(s);
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
