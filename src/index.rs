use std::collections::HashMap;

use crate::Kmer;

/// A simple position index for k-mers.
pub struct SimplePosIndex {
    k: usize,
    empty: Vec<usize>,
    index: HashMap<Kmer, Vec<usize>>,
}

impl SimplePosIndex {
    /// Create a new simple position index
    pub fn new(k: usize) -> SimplePosIndex {
        SimplePosIndex {
            k,
            empty: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// Add the forward strand of `seq` to the index.
    pub fn add_seq<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many_both_pos(self.k, seq, |pos, x, _y| {
            self.index.entry(x.clone()).or_insert(Vec::new()).push(pos);
        });
    }

    /// Add both strands of `seq` to the index.
    pub fn add_seq_both<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many_both_pos(self.k, seq, |pos, x, y| {
            self.index.entry(x.clone()).or_insert(Vec::new()).push(pos);
            self.index.entry(y.clone()).or_insert(Vec::new()).push(pos);
        });
    }

    /// Return the locations where `x` occurs.
    pub fn find(&self, x: &Kmer) -> &[usize] {
        match self.index.get(x) {
            None => &self.empty,
            Some(positions) => &positions,
        }
    }

    /// Return the possible locations for `seq` with the number of k-mers supporting each.
    pub fn find_locations<S>(&self, seq: &S) -> HashMap<i64, usize>
    where
        S: AsRef<[u8]>,
    {
        let mut res = HashMap::new();
        Kmer::with_many_both_pos(self.k, seq, |pos, x, _y| {
            let p = pos as i64;
            for q in self.find(x) {
                let r = (*q as i64) - p;
                *res.entry(r).or_insert(0) += 1;
            }
        });
        res
    }
}

/// A simple index allowing multiple sequences to be distinguished in the index.
pub struct SimpleSeqPosIndex {
    k: usize,
    name_vector: Vec<String>,
    index: HashMap<Kmer,Vec<(usize,usize)>>,
    empty: Vec<(usize, usize)>
}

impl SimpleSeqPosIndex {
    /// Create a simple index storing positions for multiple sequences.
    pub fn new(k: usize) -> SimpleSeqPosIndex {
        SimpleSeqPosIndex { k, name_vector: Vec::new(), index: HashMap::new(), empty: Vec::new() }
    }

    /// Return the vector of sequence names
    pub fn names(&self) -> &[String] {
        &self.name_vector
    }

    /// Add the forward strand of `seq` to the index.
    pub fn add_seq<S>(&mut self, name: &str, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        let n = self.name_vector.len();
        self.name_vector.push(name.to_string());

        Kmer::with_many_both_pos(self.k, seq, |pos, x, _y| {
            self.index.entry(x.clone()).or_insert(Vec::new()).push((n,pos));
        });
    }

    /// Add both strands of `seq` to the index.
    pub fn add_seq_both<S>(&mut self, name: &str, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        let n = self.name_vector.len();
        self.name_vector.push(name.to_string());

        Kmer::with_many_both_pos(self.k, seq, |pos, x, y| {
            self.index.entry(x.clone()).or_insert(Vec::new()).push((n, pos));
            self.index.entry(y.clone()).or_insert(Vec::new()).push((n, pos));
        });
    }

    /// Return the locations where `x` occurs.
    pub fn find(&self, x: &Kmer) -> &[(usize,usize)] {
        match self.index.get(x) {
            None => &self.empty,
            Some(positions) => &positions,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1() {
        let k = 11;
        let mut idx = SimplePosIndex::new(k);
        {
            let seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
            idx.add_seq(&seq);
        }
        let seq = "CACCTTAAAAACAGAGAAAGAAAAAATGACACATTACCTATACTGCAAAAACAA";
        let hits = idx.find_locations(&seq);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits.get(&46), Some(&23));
        assert_eq!(hits.get(&47), Some(&15));
    }

    #[test]
    fn test_2() {
        let k = 11;
        let mut idx = SimplePosIndex::new(k);
        {
            let seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
            idx.add_seq_both(&seq);
        }
        let seq = "CACCTTAAAAACAGAGAAAGAAAAAATGACACATTACCTATACTGCAAAAACAA";
        let hits = idx.find_locations(&seq);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits.get(&46), Some(&23));
        assert_eq!(hits.get(&47), Some(&15));
    }
}