use pyo3::prelude::*;
use std::collections::HashSet;

#[pyfunction]
fn sum_as_string(a: i32, b: i32) -> PyResult<i32> {
    Ok(a + b)
}

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
fn hamming_distance(text: &str, pattern: &str) ->  PyResult<Vec<usize>>  {
    let mut distances: Vec<usize> = Vec::new();
    let positions = pattern_count_positions_rust(text,pattern);
    let mut before = positions[0];
    for pos in &positions[1..] {
        distances.push(pos-before);
        before=*pos;
    }
    Ok(distances)
}

/// A Python module implemented in Rust.
#[pymodule]
fn bioinformatics(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hamming_distance, m)?)?;
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_count, m)?)?;
    m.add_function(wrap_pyfunction!(frequent_words, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_count_positions, m)?)?;

    Ok(())
}
