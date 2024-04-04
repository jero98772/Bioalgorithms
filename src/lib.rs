mod libs2;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use std::collections::HashSet;
use std::collections::HashMap;

use rand::seq::SliceRandom;
use rand::thread_rng;
use rand::Rng;
use rand::prelude::*;

use libs2::functions::{pattern_to_number_rust,pattern_count_frequent_words,pattern_count_positions_rust,hamming_distance,approx,generate_kmer_neighbors,d,product,log_prob,most_probable_rust,calculate_profile_matrix_greddy,score_mofit_greddy,score_mofit_random,profile_mofit_random,get_motifs,random_kmer,score_gibbs,profile_gibbs,probability_kmer,generate_gibbs,drop_one_motif};
use libs2::functions::{BASES};

#[pyfunction]
fn pattern_count(text: &str, pattern: &str) -> PyResult<i32> {
    let mut count = 0;
    let pattern_size = pattern.len();
    for i in 0..text.len() {
        let mut pattern_size = pattern_size;
        for j in 0..pattern.len() {
            if text.chars().nth(i + j).unwrap() == pattern.chars().nth(j).unwrap() {
                pattern_size -= 1;
            } else {
                break;
            }
            if pattern_size == 0 {
                count += 1;
            }
        }
    }
    Ok(count)
}


#[pyfunction]
fn frequent_words(text: &str, k: usize) -> PyResult<HashSet<&str>> {
    let mut frequent_patterns = HashSet::new();
    let mut count = vec![];
    let mut maxk = 0;

    for i in 0..text.len() - k {
        let pattern = &text[i..i + k];
        let pattern_countval = pattern_count_frequent_words(text, pattern);
        if pattern_countval > maxk {
            maxk = pattern_countval;
        }
        count.push(pattern_countval);
    }

    for i in 0..text.len() - k {
        if count[i] == maxk || count[i] > 1 {
            frequent_patterns.insert(&text[i..i + k]);
        }
    }

    println!("{}", maxk);
    Ok(frequent_patterns)
}

#[pyfunction]
fn pattern_count_positions(text: &str, pattern: &str) -> PyResult<Vec<usize>>  {
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
    Ok(positions)
}
#[pyfunction]
fn hamming_distances(text: &str, pattern: &str) ->  PyResult<Vec<usize>>  {//this function can be opticed
    let mut distances: Vec<usize> = Vec::new();
    let positions = pattern_count_positions_rust(text,pattern);
    let mut before = positions[0];
    for pos in &positions[1..] {
        distances.push(pos-before);
        before=*pos;
    }
    Ok(distances)
}

#[pyfunction]
fn clump_finding(genome: &str, k: usize, l: usize, t: usize) -> PyResult<HashSet<String>> {
    let mut clumps: HashSet<String> = HashSet::new();
    
    for i in 0..genome.len() - l + 1 {
        let window = &genome[i..i+l];
        let mut kmer_counts: HashMap<String, usize> = HashMap::new();
        
        for j in 0..window.len() - k + 1 {
            let kmer = &window[j..j+k];
            *kmer_counts.entry(kmer.to_string()).or_insert(0) += 1;
        }
        
        for (kmer, count) in kmer_counts.iter() {
            if *count >= t {
                clumps.insert(kmer.clone());
            }
        }
    }
    
    Ok(clumps)
}

#[pyfunction]
fn min_skew(text: &str) -> PyResult<(i32,usize)>{
    let mut min=std::i32::MAX;
    let mut count:i32=0;
    let mut pos=0;
    let mut min_pos=0;
    let mut c;
    for i in text.chars() {
        c=i.to_lowercase().next().unwrap();
        if c=='g'{
            count+=1;
        }else if c=='c'{
            count-=1;
        }
        if count<min{
            min=count;
            min_pos=pos;
        }
        pos+=1;
    }
    return Ok((min,min_pos));   
}


#[pyfunction]
fn approximate_pattern_matching(text: &str, pattern: &str, d: usize) -> PyResult<Vec<usize>> {
    let mut starting_positions = Vec::new();

    for i in 0..=text.len() - pattern.len() {
        if hamming_distance(&text[i..i + pattern.len()],&pattern) <= d {
            starting_positions.push(i);
        }
    }
    Ok(starting_positions)
}

#[pyfunction]
fn approximate_pattern_count(text: &str,pattern: &str, d: usize) -> PyResult<usize> {
    let mut count = 0;

    for i in 0..=text.len() - pattern.len() {
        if hamming_distance(&text[i..i + pattern.len()],&pattern) <= d {
            count += 1;
        }
    }
    Ok(count)
}

#[pyfunction]
fn frequent_words_mismatch(dna: &str, k: usize, n: usize) -> PyResult<String> {
    let mut counts = HashMap::new();
    for i in 0..=(dna.len() - k) {
        let kmer = &dna[i..(i + k)];
        *counts.entry(kmer.to_string()).or_insert(0) += 1;
    }

    let mut update_counts = HashMap::new();
    for (a, _) in &counts {
        let mut c = 0;
        for (b, _) in &counts {
            if approx(a, b, k, n) {
                c += counts.get(b).unwrap_or(&0);
            }
        }
        update_counts.insert(a.clone(), c);
    }

    let frequent = *update_counts.values().max().unwrap_or(&0);
    let mut ans = Vec::new();
    for (k, v) in &update_counts {
        if *v == frequent {
            ans.push(k.clone());
        }
    }

    Ok(ans.join(" "))
}

#[pyfunction]
fn reverse_complement(pattern: &str) -> PyResult<String> {
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
    Ok(new_chain.iter().collect::<String>())
}

#[pyfunction]
fn generate_frequency_array(text: &str, k: usize) -> PyResult<Vec<usize>> {
    let mut frequencies = vec![0; 4_usize.pow(k as u32)];
    for i in 0..=text.len() - k {
        frequencies[pattern_to_number_rust(&text[i..i + k])] += 1;
    }
    Ok(frequencies)
}
#[pyfunction]
fn pattern_to_number(kmer: &str) -> PyResult<usize> {
    let mut n = 0;
    for letter in kmer.chars() {
        n *= 4;
        n += BASES.find(letter).unwrap();
    }
    Ok(n)
}
#[pyfunction]
fn number_to_pattern(mut n: usize, k: usize) -> PyResult<String> {
    let mut pattern = String::new();
    for _ in 0..k {
        pattern.push(BASES.chars().nth(n % 4).unwrap());
        n /= 4;
    }
    Ok(pattern.chars().rev().collect::<String>())
}

#[pyfunction]
fn enumerate_motifs(_py: Python, dna: Vec<&str>, k: usize, d: usize) -> PyResult<Vec<String>> {
    let mut patterns = HashSet::new();
    for dna_string in &dna {
        for i in 0..=dna_string.len() - k {
            let kmer = &dna_string[i..i + k];
            let neighbors = generate_kmer_neighbors(kmer, d);
            for neighbor in &neighbors {
                let mut found_in_all = true;
                for dna_string2 in &dna {
                    if (0..=dna_string2.len() - k).all(|j| {
                        (0..k).all(|_l| hamming_distance(&neighbor, &dna_string2[j..j + k]) > d)
                    }) {
                        found_in_all = false;
                        break;
                    }
                }
                if found_in_all {
                    patterns.insert(neighbor.clone());
                }
            }
        }
    }
    Ok(patterns.into_iter().collect())
}

#[pyfunction]
fn median_string(_py: Python, dna: Vec<&str>, k: usize) -> PyResult<String> {
    let mut distance = usize::MAX;
    let mut median = String::new();

    for pattern in product(&['A', 'C', 'G', 'T'], k) {
        let pattern: String = pattern.iter().collect();
        let pattern = &pattern;

        if distance > d(pattern, &dna) {
            distance = d(pattern, &dna);
            median = pattern.to_string();
        }
    }

    Ok(median)
}


#[pyfunction]
fn most_probable(py: Python, text: &str, n: usize, profile: Vec<Vec<f64>>) -> PyResult<String> {
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

    Ok(most_probable_kmer)
}
#[pyfunction]
fn greedy_motif_search(py: Python, k: usize, t: usize, dna: Vec<String>, pseudo_counts: Option<bool>) -> PyResult<Vec<String>> {
    let pseudo_counts = match pseudo_counts {
        Some(value) => value,
        None => false,
    };

    let mut best_motifs = dna.iter().map(|genome| genome.chars().take(k).collect()).collect::<Vec<String>>();
    for motif in dna[0].chars().take(dna[0].len() - k + 1) {
        let mut motifs = vec![motif.to_string()];
        for i in 1..t {
            motifs.push(most_probable_rust(&dna[i], k, calculate_profile_matrix_greddy(&motifs, BASES, k, pseudo_counts)));
        }
        if score_mofit_greddy(&motifs, BASES, k, pseudo_counts) < score_mofit_greddy(&best_motifs, BASES, k, pseudo_counts) {
            best_motifs = motifs;
        }
    }

    Ok(best_motifs)
}

#[pyfunction]
fn randomized_motif_search(py: Python, k: usize, t: usize, dna: Vec<String>, eps: i32) -> PyResult<(i32, Vec<String>)> {
    let bases = "ACGT";
    let mut best_motifs = vec![];

    for i in 0..t {
        let mut rng = thread_rng();
        let mut motifs: Vec<String> = dna.iter().map(|s| random_kmer(&mut rng, s, k)).collect();
        best_motifs = motifs.clone();

        loop {
            let profile = profile_mofit_random(&motifs, k, eps);
            motifs = get_motifs(&profile, &dna, k);
            if score_mofit_random(k, &motifs, &bases) < score_mofit_random(k, &best_motifs, &bases) {
                best_motifs = motifs.clone();
            } else {
                return Ok((score_mofit_random(k, &best_motifs, &bases), best_motifs));
            }
        }
    }

    Ok((score_mofit_random(k, &best_motifs, &bases), best_motifs))
}

#[pyfunction]
fn randomized_motif_search_driver(py: Python, k: usize, t: usize, dna: Vec<String>, n: usize) -> PyResult<(f64, Vec<String>)> {
    let mut best = std::f64::MAX;
    let mut mm = vec![];

    for i in 0..n {
        let (sc, motifs) = randomized_motif_search(py, k, t, dna.clone(),n as i32)?;
        let sc_f64 = sc as f64; // Convert score to f64
        if sc_f64 < best {
            best = sc_f64;
            mm = motifs.clone();
        }
        if i % 100 == 0 {
            println!("{}, {}", i, best);
            for m in &mm {
                println!("{}", m);
            }
        }
    }
    Ok((best, mm))
}
#[pyfunction]
fn gibbs(py: Python, k: usize, t: usize, n: usize, dna: Vec<String>, eps: i32) -> PyResult<(i32, Vec<String>, Vec<i32>)> {
    let mut best_score = std::i32::MAX;
    let mut best_motifs = vec![];
    let mut rng = rand::thread_rng();
    let mut trace = vec![];
    let mut motifs = Vec::new();
    for i in 0..t {
        motifs.push(random_kmer(&mut rng, &dna[i], k));
    }
    best_motifs = motifs.clone();

    for _ in 0..n {
        let i =  rng.gen_range(0..t);
        let profile = profile_gibbs(&drop_one_motif(&motifs, i), BASES, k, eps);
        let probabilities: Vec<f64> = (0..(dna[i].len() - k + 1)).map(|ll| probability_kmer(&dna[i][ll..(ll + k)], &profile, BASES)).collect();
        let motif_index = generate_gibbs(&probabilities);
        motifs[i] = dna[i][motif_index..(motif_index + k)].to_string();
        let sc = score_gibbs(k, &motifs, BASES);
        if sc < best_score {
            best_score = sc;
            best_motifs = motifs.clone();
        }
        trace.push(best_score);
    }

    Ok((best_score, best_motifs, trace))
}

#[pyfunction]
fn distance_between_pattern_and_strings(_py: Python, pattern: &str, dna: Vec<&str>) -> usize {
    dna.iter().map(|motif| {
        (0..(motif.len() - pattern.len() + 1)).map(|i| {
            hamming_distance(pattern, &motif[i..(i + pattern.len())])
        }).min().unwrap()
    }).sum()
}

#[pyfunction]
fn kmer_composition(_py: Python, k: usize, dna: &str) -> Vec<String> {
    (0..=dna.len() - k).map(|i| dna[i..i+k].to_string()).collect()
}

#[pymodule]
fn bioinformatics(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(kmer_composition, m)?)?;
    m.add_function(wrap_pyfunction!(distance_between_pattern_and_strings, m)?)?;
    m.add_function(wrap_pyfunction!(gibbs, m)?)?;
    m.add_function(wrap_pyfunction!(randomized_motif_search, m)?)?;
    m.add_function(wrap_pyfunction!(greedy_motif_search, m)?)?;
    m.add_function(wrap_pyfunction!(most_probable, m)?)?;
    m.add_function(wrap_pyfunction!(median_string, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_motifs, m)?)?;
    m.add_function(wrap_pyfunction!(number_to_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_to_number, m)?)?;
    m.add_function(wrap_pyfunction!(generate_frequency_array, m)?)?;
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(frequent_words_mismatch, m)?)?;
    m.add_function(wrap_pyfunction!(approximate_pattern_matching, m)?)?;
    m.add_function(wrap_pyfunction!(approximate_pattern_count, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_distances, m)?)?;
    m.add_function(wrap_pyfunction!(min_skew, m)?)?;
    m.add_function(wrap_pyfunction!(clump_finding, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_count, m)?)?;
    m.add_function(wrap_pyfunction!(frequent_words, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_count_positions, m)?)?;

    Ok(())
}
