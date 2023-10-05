use std::collections::BTreeMap;

use crate::{
    counting_kmer_frequency_iterator, frequency::frequency_vector_iter, merge, unique, Kmer,
};

/// Accumulate k-mer frequencies for small to medium k.
pub struct DenseAccumulator {
    k: usize,
    counts: Vec<usize>,
}

impl DenseAccumulator {
    /// Construct a new accumulator.
    pub fn new(k: usize) -> DenseAccumulator {
        let n = 1usize << (2 * k);
        let mut counts = Vec::new();
        counts.resize(n, 0);
        DenseAccumulator { k, counts }
    }

    /// Add k-mers from the forward strand of the given sequence.
    pub fn add<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many(self.k, seq, |x| {
            self.counts[x.0 as usize] += 1;
        });
    }

    /// Add k-mers from both strands of the given sequence.
    pub fn add_both<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many_both(self.k, seq, |x, y| {
            self.counts[x.0 as usize] += 1;
            self.counts[y.0 as usize] += 1;
        });
    }

    /// Get a k-mer frequency iterator over the accumulator.
    pub fn kmer_frequencies(&self) -> impl Iterator<Item = (Kmer, usize)> + '_ {
        frequency_vector_iter(self.counts.iter())
    }
}

/// A simple sparse accumulator
pub struct SimpleSparseAccumulator {
    k: usize,
    counts: BTreeMap<Kmer, usize>,
}

impl SimpleSparseAccumulator {
    /// Construct a new sparse kmer frequency accumulator
    pub fn new(k: usize) -> SimpleSparseAccumulator {
        SimpleSparseAccumulator {
            k,
            counts: BTreeMap::new(),
        }
    }

    /// Add k-mers from the forward strand of the given sequence.
    pub fn add<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many(self.k, seq, |x| {
            *self.counts.entry(x.clone()).or_insert(0) += 1;
        });
    }

    /// Add k-mers from both strands of the given sequence.
    pub fn add_both<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many_both(self.k, seq, |x, y| {
            *self.counts.entry(x.clone()).or_insert(0) += 1;
            *self.counts.entry(y.clone()).or_insert(0) += 1;
        });
    }

    /// Get a k-mer frequency iterator over the accumulator.
    pub fn kmer_frequencies(&self) -> impl Iterator<Item = (Kmer, usize)> + '_ {
        self.counts.iter().map(|(x, f)| (x.clone(), *f))
    }
}

/// A sparse k-mer accumulator for large set of k-mers.
pub struct BigSparseAccumulator {
    k: usize,
    buffer_size: usize,
    buffer: Vec<Kmer>,
    unique_kmers: Vec<Kmer>,
    non_unique_kmers: Vec<(Kmer, usize)>,
}

impl BigSparseAccumulator {
    /// Create a new sparse k-mer accumulator, buffering `buffer_size` k-mers between compactions.
    pub fn new(k: usize, buffer_size: usize) -> BigSparseAccumulator {
        BigSparseAccumulator {
            k,
            buffer_size,
            buffer: Vec::with_capacity(buffer_size),
            unique_kmers: Vec::new(),
            non_unique_kmers: Vec::new(),
        }
    }

    /// Add k-mers from the forward strand of the given sequence.
    pub fn add<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many(self.k, seq, |x| {
            self.buffer.push(x.clone());
            if self.buffer.len() + 2 >= self.buffer_size {
                self.flush();
            }
        });
    }

    /// Add k-mers from both strands of the given sequence.
    pub fn add_both<S>(&mut self, seq: &S)
    where
        S: AsRef<[u8]>,
    {
        Kmer::with_many_both(self.k, seq, |x, y| {
            self.buffer.push(x.clone());
            self.buffer.push(y.clone());
            if self.buffer.len() + 2 >= self.buffer_size {
                self.flush();
            }
        });
    }

    /// Get a k-mer frequency iterator over the accumulator.
    pub fn kmer_frequencies(&mut self) -> impl Iterator<Item = (Kmer, usize)> + '_ {
        self.flush();
        let unique_kmer_frequencies = unique(self.unique_kmers.iter());
        let non_unique_kmer_frequencies =
            self.non_unique_kmers.iter().map(|(x, f)| (x.clone(), *f));
        merge(unique_kmer_frequencies, non_unique_kmer_frequencies)
    }

    fn flush(&mut self) {
        if self.buffer.len() == 0 {
            return;
        }

        self.buffer.sort();

        let mut unique_kmers = Vec::new();
        let mut non_unique_kmers = Vec::new();
        {
            let new_kmer_frequencies = counting_kmer_frequency_iterator(self.buffer.iter());
            let unique_kmer_frequencies = unique(self.unique_kmers.iter());
            let non_unique_kmer_frequencies =
                self.non_unique_kmers.iter().map(|(x, f)| (x.clone(), *f));

            for (x, f) in merge(
                merge(new_kmer_frequencies, unique_kmer_frequencies),
                non_unique_kmer_frequencies,
            ) {
                if f == 1 {
                    unique_kmers.push(x);
                } else {
                    non_unique_kmers.push((x, f));
                }
            }
        }
        self.buffer.clear();
        std::mem::swap(&mut unique_kmers, &mut self.unique_kmers);
        std::mem::swap(&mut non_unique_kmers, &mut self.non_unique_kmers);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MiniRng {
        x: u64,
    }

    impl MiniRng {
        fn new(seed: u64) -> MiniRng {
            MiniRng { x: seed }
        }

        fn rnd(&mut self) -> u64 {
            self.x = self.x.wrapping_mul(2862933555777941757u64);
            self.x = self.x.wrapping_add(3037000493u64);
            self.x
        }
    }

    fn make_some_sequence(rng: &mut MiniRng) -> String {
        let mut s = String::new();
        let k: usize = 17;
        let m = (1u64 << (2*k)) - 1;
        for _i in 0..20 {
            let u = rng.rnd() & m;
            let x = Kmer::from_u64(u);
            s += &x.render(k);
        }
        println!("{}", s);
        s
    }

    #[test]
    fn test_1() {
        let k = 7;

        let seq1 = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
        let seq2 = "TTATATGCACATGATTTACCTCCAGTAAAGATTTCTAAAGTAAATGTAAGAAAAGGGACAAGTTCTCCAGGGAAAGAAAACCTGACAACAAAGTGTGTCTAAACTTGTGATCAAAGGTTAAATATCTACATCAAAGCCCTAATCAGCACA";

        let counts1 = Kmer::frequency_vector(k, &seq1);
        let counts2 = Kmer::frequency_vector_both(k, &seq2);
        assert_eq!(counts1.len(), counts2.len());

        let mut acc = DenseAccumulator::new(k);
        acc.add(&seq1);
        {
            let mut i = 0;
            for (x, f) in acc.kmer_frequencies() {
                assert!(x.0 < counts1.len() as u64);
                while (i as u64) < x.0 {
                    assert_eq!(counts1[i], 0);
                    i += 1;
                }
                assert_eq!(f, counts1[i]);
                i += 1;
            }
            assert!(i <= counts1.len());
        }
        acc.add_both(&seq2);
        {
            let mut i = 0;
            for (x, f) in acc.kmer_frequencies() {
                assert!(x.0 < counts1.len() as u64);
                while (i as u64) < x.0 {
                    assert_eq!(counts1[i] + counts2[i], 0);
                    i += 1;
                }
                assert_eq!(f, counts1[i] + counts2[i]);
                i += 1;
            }
            assert!(i <= counts1.len());
        }
    }

    #[test]
    fn test_2() {
        let k = 7;

        let seq1 = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
        let seq2 = "TTATATGCACATGATTTACCTCCAGTAAAGATTTCTAAAGTAAATGTAAGAAAAGGGACAAGTTCTCCAGGGAAAGAAAACCTGACAACAAAGTGTGTCTAAACTTGTGATCAAAGGTTAAATATCTACATCAAAGCCCTAATCAGCACA";

        let counts1 = Kmer::frequency_vector(k, &seq1);
        let counts2 = Kmer::frequency_vector_both(k, &seq2);
        assert_eq!(counts1.len(), counts2.len());

        let mut acc = SimpleSparseAccumulator::new(k);
        acc.add(&seq1);
        {
            let mut i = 0;
            for (x, f) in acc.kmer_frequencies() {
                assert!(x.0 < counts1.len() as u64);
                while (i as u64) < x.0 {
                    assert_eq!(counts1[i], 0);
                    i += 1;
                }
                assert_eq!(f, counts1[i]);
                i += 1;
            }
            assert!(i <= counts1.len());
        }
        acc.add_both(&seq2);
        {
            let mut i = 0;
            for (x, f) in acc.kmer_frequencies() {
                assert!(x.0 < counts1.len() as u64);
                while (i as u64) < x.0 {
                    assert_eq!(counts1[i] + counts2[i], 0);
                    i += 1;
                }
                assert_eq!(f, counts1[i] + counts2[i]);
                i += 1;
            }
            assert!(i <= counts1.len());
        }
    }

    #[test]
    fn test_3() {
        let k = 7;
        let mut acc = BigSparseAccumulator::new(k, 1024);
        let mut rng = MiniRng::new(19);
        for _i in 0..10 {
            let seq = make_some_sequence(&mut rng);
            acc.add(&seq.as_bytes().iter());
        }
        for _i in 0..10 {
            let seq = make_some_sequence(&mut rng);
            acc.add_both(&seq.as_bytes().iter());
        }
        acc.flush();
        let xfs = Vec::from_iter(acc.kmer_frequencies());
        assert_eq!(xfs.len(), 7585);
        let mut h = Vec::new();
        for (_x, f) in xfs.iter() {
            let c = *f;
            while h.len() <= c {
                h.push(0);
            }
            h[c] += 1;
        }
        assert_eq!(h, vec![0, 5553, 1680, 308, 37, 7]);
    }
}
