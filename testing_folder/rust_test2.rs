
fn most_probable(text: &str, n: usize, profile: &Vec<Vec<f64>>) -> String {
    let bases = "ACGT";

    fn log_prob(kmer: &str, n: usize, profile: &Vec<Vec<f64>>, bases: &str) -> f64 {
        kmer.chars().enumerate().map(|(j, c)| profile[bases.find(c).unwrap()][j]).fold(0.0, |acc, p| acc + p.ln())
    }

    let find_most_probable = |text: &str, n: usize, profile: &Vec<Vec<f64>>, bases: &str| {
        let mut kmers = Vec::new();
        for i in 0..=text.len() - n {
            kmers.push(&text[i..i + n]);
        }
        kmers.iter().max_by(|&&a, &&b| {
            log_prob(a, n, profile, bases)
                .partial_cmp(&log_prob(b, n, profile, bases))
                .unwrap()
        }).unwrap().to_string()
    };

    find_most_probable(text, n, profile, bases)
}

fn main() {
    // Example usage
    let text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT";
    let n = 5;
    let profile = vec![
        vec![0.2, 0.2, 0.3, 0.2, 0.3],
        vec![0.4, 0.3, 0.1, 0.5, 0.1],
        vec![0.3, 0.3, 0.5, 0.2, 0.4],
        vec![0.1, 0.2, 0.1, 0.1, 0.2],
    ];
    let result = most_probable(text, n, &profile);
    println!("{}", result); // Output: CCGAG
}

