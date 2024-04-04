
pub mod functions {
    use rand::Rng;
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

    pub fn calculate_profile_matrix_greddy(motifs: &[String], bases: &str, k: usize, pseudo_counts: bool) -> Vec<Vec<f64>> {
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

    pub fn score_mofit_random(k: usize, motifs: &[String], bases: &str) -> i32 {
        let mut total = 0;
        for j in 0..k {
            let mut counts = vec![0; bases.len()];
            for motif in motifs {
                let i = bases.find(motif.chars().nth(j).unwrap()).unwrap();
                counts[i] += 1;
            }
            let mut max = -1;
            let mut ii = -1;
            for i in 0..bases.len() {
                if max < counts[i] {
                    ii = i as i32;
                    max = counts[ii as usize];
                }
            }
            for i in 0..bases.len() {
                if (i as i32) != ii {
                    total += counts[i];
                }
            }
        }
        total
    }
    pub fn profile_mofit_random(motifs: &[String], k: usize, eps: i32) -> Vec<Vec<f64>> {
        let mut matrix = vec![vec![eps as f64; k]; 4];
        for kmer in motifs {
            for j in 0..k {
                let i = "ACGT".find(kmer.chars().nth(j).unwrap()).unwrap();
                matrix[i][j] += 1.0;
            }
        }
        let motifs_len = motifs.len() as f64;
        matrix.iter_mut().for_each(|row| row.iter_mut().for_each(|x| *x /= motifs_len));
        matrix
    }
    pub fn get_motifs(profile: &[Vec<f64>], dna: &[String], k: usize) -> Vec<String> {
        let mut motifs = vec![];
        for s in dna {
            let mut max_probability = -1.0;
            let mut most_probable_kmer = String::new();
            for kmer in s.chars().collect::<Vec<char>>().windows(k) {
                let kmer_str: String = kmer.iter().collect();
                let prob = kmer.iter().enumerate().fold(1.0, |p, (j, ch)| {
                    let i = "ACGT".find(*ch).unwrap();
                    p * profile[i][j]
                });
                if prob > max_probability {
                    max_probability = prob;
                    most_probable_kmer = kmer_str;
                }
            }
            motifs.push(most_probable_kmer);
        }
        motifs
    }

    pub fn random_kmer(rng: &mut rand::rngs::ThreadRng, string: &str, k: usize) -> String {
        let i = rng.gen_range(0..string.len() - k + 1);
        string[i..i + k].to_string()
    }
    pub fn score_gibbs(k: usize, motifs: &[String], bases: &str) -> i32 {
        let mut total = 0;
        for j in 0..k {
            let mut counts = vec![0; bases.len()];
            for motif in motifs {
                let i = bases.find(motif.chars().nth(j).unwrap()).unwrap();
                counts[i] += 1;
            }
            let mut max = -1;
            let mut ii = -1;
            for i in 0..bases.len() {
                if max < counts[i] {
                    ii = i as i32;
                    max = counts[ii as usize];
                }
            }
            for i in 0..bases.len() {
                if (i as i32) != ii {
                    total += counts[i];
                }
            }
        }
        total
    }
    pub fn counts(motifs: &[String], bases: &str, k: usize, eps: i32) -> Vec<Vec<i32>> {
        let mut matrix = vec![vec![eps; k]; bases.len()];
        for kmer in motifs {
            for (j, ch) in kmer.chars().enumerate() {
                let i = bases.find(ch).unwrap();
                matrix[i][j] += 1;
            }
        }
        matrix
    }
    pub fn profile_gibbs(motifs: &[String], bases: &str, k: usize, eps: i32) -> Vec<Vec<f64>> {
        let matrix = counts(motifs, bases, k, eps);
        let motifs_len = motifs.len() as f64;
        matrix.iter().map(|row| {
            row.iter().map(|&count| count as f64 / motifs_len).collect()
        }).collect()
    }
    pub fn probability_kmer(kmer: &str, profile: &[Vec<f64>], bases: &str) -> f64 {
        let mut p = 1.0;
        for (j, ch) in kmer.chars().enumerate() {
            let i = bases.find(ch).unwrap();
            p *= profile[i][j];
        }
        p
    }
    pub fn generate_gibbs(probabilities: &[f64]) -> usize {
        let rr = rand::thread_rng().gen::<f64>();
        let mut accumulated = 0.0;
        for (i, &p) in probabilities.iter().enumerate() {
            accumulated += p;
            if rr <= accumulated {
                return i;
            }
        }
        probabilities.len() - 1
    }
    pub fn drop_one_motif(motifs: &[String], i: usize) -> Vec<String> {
        motifs.iter().enumerate().filter_map(|(j, motif)| if j != i { Some(motif.clone()) } else { None }).collect()
    }

}
