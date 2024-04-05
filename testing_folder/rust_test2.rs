
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

use std::collections::HashMap;

fn kmer_composition(k: usize, text: &str) -> Vec<&str> {
    text.chars().collect::<Vec<char>>()
        .windows(k)
        .map(|window| window.iter().collect::<String>())
        .collect::<Vec<String>>()
}

fn grph_kmers(strings: &[&str]) -> Vec<(&str, &str)> {
    let kk = strings[0].len() - 1;
    let mut graph = Vec::new();

    for &s in strings {
        for &t in strings {
            if s != t && &s[s.len() - kk..] == &t[0..kk] {
                graph.push((s, t));
            }
        }
    }

    graph
}

fn de_bruijn(k: usize, text: &str) -> Vec<(&str, Vec<&str>)> {
    let kmers = kmer_composition(k - 1, text);

    fn standardize(ll: Vec<&str>) -> Vec<&str> {
        let mut lll = ll.to_vec();
        lll.sort();
        lll.dedup();
        lll
    }

    let pathgraph = grph_kmers(&kmers);

    let mut de_bruijn_dict: HashMap<&str, Vec<&str>> = HashMap::new();
    for (a, b) in pathgraph {
        de_bruijn_dict.entry(a).or_insert(Vec::new()).push(b);
    }

    let mut graph: Vec<(&str, Vec<&str>)> = de_bruijn_dict.iter()
        .map(|(a, b)| (*a, standardize(b.to_vec())))
        .collect();
    graph.sort();

    graph
}



fn main() {
    // Example usage
    let k = 4;
    let text = "AAGATTCTCTAAGA";

    let de_bruijn_graph = de_bruijn(k, text);
    println!("{:?}", de_bruijn_graph);


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

