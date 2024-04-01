
pub mod functions {
    use std::collections::HashSet;

    pub const BASES: &str = "ACGT";
    pub fn pattern_to_number_rust(kmer: &str) -> usize {
        let mut n = 0;
        for letter in kmer.chars() {
            n *= 4;
            n += BASES.find(letter).unwrap();
        }
        n
    }
    

    pub fn pattern_count_frequent_words(text: &str, pattern: &str) -> usize {
        let mut count = 0;
        let pattern_size = pattern.len();

        for i in 0..text.len() {
            let mut pattern_size = pattern_size;
            for j in 0..pattern.len() {
                if j < pattern.len() && i + j < text.len() && text.chars().nth(i + j).unwrap() == pattern.chars().nth(j).unwrap() {
                    pattern_size -= 1;
                } else {
                    break;
                }
                if pattern_size == 0 {
                    count += 1;
                }
            }
        }
        count
    }

    pub fn pattern_count_positions_rust(text: &str, pattern: &str) -> Vec<usize>  {
        //let mut count = 0;
        let mut positions: Vec<usize> = Vec::new();
        let pattern_size = pattern.len();

        for i in 0..text.len() {
            let mut pattern_size = pattern_size;
            for j in 0..pattern.len() {
                if j < pattern.len() && i + j < text.len() && text.chars().nth(i + j).unwrap() == pattern.chars().nth(j).unwrap() {
                    pattern_size -= 1;
                } else {
                    break;
                }
                if pattern_size == 0 {
                    positions.push(i);
                    //count += 1;
                }
            }
        }
        positions
    }
    pub fn hamming_distance(str1: &str, str2: &str) -> usize {
        let mut mismatch = 0;

        for (char1, char2) in str2.chars().zip(str1.chars()) {
            if char1 != char2 {
                mismatch += 1;
            }
        }
        mismatch
    }


    pub fn approx(a: &str, b: &str, _k: usize, n: usize) -> bool {
        let mut mismatch = 0;
        for (ca, cb) in a.chars().zip(b.chars()) {
            if ca != cb {
                mismatch += 1;
            }
            if mismatch > n {
                return false;
            }
        }
        true
    }

    pub fn reverse_pattern(pattern: &str) -> String {
            let chain = pattern.to_lowercase();
            //let chain = "ATGATCAAG";
            let mut new_chain = Vec::new();
            for char in chain.chars().rev() {
                match char.to_ascii_lowercase() {
                    'a' => new_chain.push('t'),
                    't' => new_chain.push('a'),
                    'g' => new_chain.push('c'),
                    'c' => new_chain.push('g'),
                    _ => {}
                }
            }
            new_chain.iter().collect::<String>()
        }
    pub fn generate_kmer_neighbors(pattern: &str, d: usize) -> HashSet<String> {
        let alphabet = ['A', 'C', 'G', 'T'];
        let mut neighbors = HashSet::new();
        generate_neighbors_rust(pattern.as_bytes(), d, &alphabet, &mut neighbors);
        neighbors
    }

    pub fn generate_neighbors_rust(pattern: &[u8], d: usize, alphabet: &[char], neighbors: &mut HashSet<String>) {
        if d == 0 {
            neighbors.insert(String::from_utf8_lossy(pattern).to_string());
            return;
        }
        for i in 0..pattern.len() {
            let original = pattern[i];
            for &ch in alphabet {
                if ch as u8 != original {
                    let mut mutated_pattern = pattern.to_vec();
                    mutated_pattern[i] = ch as u8;
                    generate_neighbors_rust(&mutated_pattern, d - 1, alphabet, neighbors);
                }
            }
        }
    }
    pub fn d(pattern: &str, dna: &[&str]) -> usize {
        dna.iter()
            .map(|dna_seq| {
                (0..dna_seq.len() - pattern.len() + 1)
                    .map(|i| hamming_distance(pattern, &dna_seq[i..i + pattern.len()]))
                    .min()
                    .unwrap_or(0)
            })
            .sum()
    }

    pub fn product<T: Clone>(choices: &[T], repeat: usize) -> Vec<Vec<T>> {
        if repeat == 0 {
            vec![vec![]]
        } else {
            let base = product(choices, repeat - 1);
            let mut result = Vec::new();
            for c in choices {
                for p in base.iter() {
                    let mut v = p.clone();
                    v.push(c.clone());
                    result.push(v);
                }
            }
            result
        }
    }
    pub fn log_prob(kmer: &str, profile: &[Vec<f64>]) -> f64 {
        let mut sum = 0.0;
        for (i, base) in kmer.chars().enumerate() {
            sum += profile[BASES.find(base).unwrap()][i].ln();
        }
        sum
    }
    pub fn most_probable_rust( text: &str, n: usize, profile: Vec<Vec<f64>>) -> String {
        let mut max_prob = f64::NEG_INFINITY;
        let mut most_probable_kmer = String::new();
        for i in 0..=(text.len() - n) {
            let kmer = &text[i..(i + n)];
            let prob = log_prob(kmer, &profile);
            if prob > max_prob {
                max_prob = prob;
                most_probable_kmer = kmer.to_string();
            }
        }

        most_probable_kmer
    }
    pub fn count_occurrences_of_bases(motifs: &[String], bases: &str, k: usize, pseudo_counts: bool) -> Vec<Vec<i32>> {
        let mut matrix = vec![vec![1; k]; bases.len()];
        if !pseudo_counts {
            matrix.iter_mut().for_each(|row| {
                row.iter_mut().for_each(|elem| *elem = 0);
            });
        }

        for kmer in motifs {
            for (j, &base) in kmer.as_bytes().iter().enumerate() {
                let i = bases.bytes().position(|x| x == base).unwrap();
                matrix[i][j] += 1;
            }
        }
        matrix
    }

    pub fn calculate_profile_matrix(motifs: &[String], bases: &str, k: usize, pseudo_counts: bool) -> Vec<Vec<f64>> {
        let matrix = count_occurrences_of_bases(motifs, bases, k, pseudo_counts);
        let total = motifs.len() as f64;
        matrix.iter().map(|row| {
            row.iter().map(|&count| count as f64 / total).collect()
        }).collect()
    }

    pub fn score_mofit_greddy(motifs: &[String], bases: &str, k: usize, pseudo_counts: bool) -> usize {
        let matrix = count_occurrences_of_bases(motifs, bases, k, pseudo_counts);
        let mut total = 0;
        for j in 0..k {
            let mut m = 0;
            for i in 0..bases.len() {
                if m < matrix[i][j] {
                    m = matrix[i][j];
                }
            }
            total += bases.len() - m as usize;
        }
        total
    }
}
