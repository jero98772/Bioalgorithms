
pub mod functions {

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


    pub fn approx(a: &str, b: &str, k: usize, n: usize) -> bool {
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

}