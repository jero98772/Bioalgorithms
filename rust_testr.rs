use std::collections::HashSet;
use std::collections::HashMap;

fn pattern_count(text: &str, pattern: &str) -> usize {
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
fn pattern_count_positions(text: &str, pattern: &str) -> Vec<usize>  {
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

fn frequent_words(text: &str, k: usize) -> HashSet<&str> {
    let mut frequent_patterns = HashSet::new();
    let mut count = vec![];
    let mut maxk = 0;

    for i in 0..text.len() - k {
        let pattern = &text[i..i + k];
        let pattern_countval = pattern_count(text, pattern);
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
    frequent_patterns
}
fn reverse_complement(text: &str)-> String{
    let mut reverse = String::new();
    let total = text.len();
    let mut c;
    let mut lowercase_c;
    for i in 0..total {
        c=text.chars().nth(total-i-1).unwrap();
        lowercase_c = c.to_lowercase().next().unwrap();
        if lowercase_c=='a'{
            reverse.push('t');
        } else if lowercase_c=='t'{
            reverse.push('a');
        } else if lowercase_c=='c'{
            reverse.push('g');
        } else if lowercase_c=='g'{
            reverse.push('c');
        }
        
    }
    reverse
}


fn hamming_distances(text: &str, pattern: &str) -> Vec<usize>  {
    let mut distances: Vec<usize> = Vec::new();
    let positions = pattern_count_positions(text,pattern);
    let mut before = positions[0];
    for pos in &positions[1..] {
        distances.push((pos-before));
        before=*pos;
    }
    distances
}
//APPROXIMATEPATTERNCOUNT

fn hamming_distance(str1: &str, str2: &str) -> usize {
    let mut mismatch = 0;

    for (char1, char2) in str2.chars().zip(str1.chars()) {
        if char1 != char2 {
            mismatch += 1;
        }
    }
    mismatch
}

fn approximate_pattern_matching(text: &str, pattern: &str, d: usize) -> Vec<usize> {
    let mut starting_positions = Vec::new();

    for i in 0..=text.len() - pattern.len() {
        if hamming_distance(&text[i..i + pattern.len()],&pattern) <= d {
            starting_positions.push(i);
        }
    }
    starting_positions
}

fn approximate_pattern_count(pattern: &str, text: &str, d: usize) -> usize {
    let mut count = 0;

    for i in 0..=text.len() - pattern.len() {
        if hamming_distance(&text[i..i + pattern.len()],&pattern) <= d {
            count += 1;
        }
    }
    count
}



fn clump_finding(genome: &str, k: usize, l: usize, t: usize) -> HashSet<String> {
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
    
    clumps
}
fn min_skew(text: &str) -> (i32,usize){
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
    return (min,min_pos);   
}


fn main() {

    let text = "AACAAGCATAAACATTAAAGAG";
    let pattern = "AAAAA";
    let d = 2;

    println!(
        "Approximate pattern count: {}",
        approximate_pattern_count(pattern, text, d)
    );

    println!(
        "Approximate pattern matching starting positions: {:?}",
        approximate_pattern_matching(text, pattern, d)
    );
    let minskew=min_skew("CATGGGCATCGGCCATACGCC");
    println!("{} {}",minskew.0,minskew.1 );
    let genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAG";
    let k = 5;
    let l = 50;
    let t = 4;
    
    let result = clump_finding(genome, k, l, t);
    println!("{:?}", result);

    let text = "AACAAGCATAAACATTAAAGAG";
    let pattern = "AAAA";
    println!("{}",approximate_pattern_count(text,pattern,2));
    let pattern_positions=hamming_distances(text,pattern);
    //println!("pattern_count_positions: {}",;
    for pos in pattern_positions {
        println!("patterns {}", pos);
    }
    println!("reverse: {}",reverse_complement(text));
    let result = frequent_words(text, 3);
    println!("Max count: {}", result.len());
    println!("Frequent patterns:");
    for pattern in result {
        println!("{}", pattern);
    }
}