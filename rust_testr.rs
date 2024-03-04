use std::collections::HashSet;

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


fn hamming_distance(text: &str, pattern: &str) -> Vec<usize>  {
    let mut distances: Vec<usize> = Vec::new();
    let positions = pattern_count_positions(text,pattern);
    let mut before = positions[0];
    for pos in &positions[1..] {
        distances.push((pos-before));
        before=*pos;
    }
    distances
}

fn main() {
    let text = "actgactcccaccccc";
    let pattern = "a";
    let pattern_positions=hamming_distance(text,pattern);
    //println!("pattern_count_positions: {}",;
    for pos in pattern_positions {
        println!("{}", pos);
    }
    println!("reverse: {}",reverse_complement(text));
    let result = frequent_words(text, 3);
    println!("Max count: {}", result.len());
    println!("Frequent patterns:");
    for pattern in result {
        println!("{}", pattern);
    }
}
