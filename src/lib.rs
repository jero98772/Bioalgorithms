use pyo3::prelude::*;
use std::collections::HashSet;
use std::collections::HashMap;


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

fn pattern_count_frequent_words(text: &str, pattern: &str) -> usize {
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
fn pattern_count_positions_rust(text: &str, pattern: &str) -> Vec<usize>  {
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

fn hamming_distance(str1: &str, str2: &str) -> usize {
    let mut mismatch = 0;

    for (char1, char2) in str2.chars().zip(str1.chars()) {
        if char1 != char2 {
            mismatch += 1;
        }
    }
    mismatch
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


/// A Python module implemented in Rust.
#[pymodule]
fn bioinformatics(_py: Python, m: &PyModule) -> PyResult<()> {
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
