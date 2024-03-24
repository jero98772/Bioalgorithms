mod libs2;

use pyo3::prelude::*;
use std::collections::HashSet;
use std::collections::HashMap;

use libs2::functions::{pattern_to_number_rust,pattern_count_frequent_words,pattern_count_positions_rust,hamming_distance,approx,generate_kmer_neighbors,d,product};
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

/// A Python module implemented in Rust.
#[pymodule]
fn bioinformatics(_py: Python, m: &PyModule) -> PyResult<()> {
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
