use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

use crate::Kmer;

/// A Trait for capturing k-mer canonicalisation.
///
/// K-mers come in pairs over reverse complementation. This trait generalises
/// an interface for implementing a variety of normalisation methods.
pub trait Canonical {
    /// Return true iff `x` is the canonical choice out of `x` and it's reverse complement.
    fn is_canonical(&self, x: &Kmer) -> bool {
        *x == self.canonical(x)
    }

    /// Return the canonical form for `x` which will be either `x` or its reverse complement.
    fn canonical(&self, x: &Kmer) -> Kmer;
}

/// A simple canonicalisation based on lexicographic ordering
pub struct CanonicalLex {
    k: usize,
}

impl CanonicalLex {
    /// Create a new lexicographic canonicalisation object.
    pub fn new(k: usize) -> CanonicalLex {
        CanonicalLex { k }
    }
}

impl Canonical for CanonicalLex {
    fn is_canonical(&self, x: &Kmer) -> bool {
        let s = 2 * (self.k - 1);
        let msb = x.0 >> s;
        let lsb = 3 - (x.0 & 3);
        println!("x = {}", x.render(self.k));
        println!("msb = {}, lsb = {}", msb, lsb);
        if msb != lsb {
            return msb < lsb;
        }
        *x == self.canonical(x)
    }

    fn canonical(&self, x: &Kmer) -> Kmer {
        let x_bar = x.rev_comp(self.k);
        if *x <= x_bar {
            x.clone()
        } else {
            x_bar
        }
    }
}

/// Canonicalisation based on hashing.
pub struct CanonicalHash {
    k: usize,
}

impl CanonicalHash {
    /// Create a new hash based canonicalisation object.
    pub fn new(k: usize) -> CanonicalHash {
        CanonicalHash { k }
    }

    /// hash a value.
    fn hash(x: &Kmer) -> u64 {
        let mut h = DefaultHasher::new();
        h.write_u64(x.0);
        h.finish()
    }
}

impl Canonical for CanonicalHash {
    fn canonical(&self, x: &Kmer) -> Kmer {
        let y = x.rev_comp(self.k);
        let x_hash = CanonicalHash::hash(x);
        let y_hash = CanonicalHash::hash(&y);
        if x_hash <= y_hash {
            x.clone()
        } else {
            y
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lex_1() {
        let k = 11;
        let seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";

        Kmer::with_many_both(k, &seq, |x, y| {
            let c = CanonicalLex::new(k);
            let x_can = c.canonical(x);
            let y_can = c.canonical(y);
            assert_eq!(x_can, y_can);
            assert_eq!(c.is_canonical(x), (*x <= *y));
            assert_eq!(c.is_canonical(y), (*y <= *x));
        });
    }

    #[test]
    fn test_hash_1() {
        let k = 11;
        let seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";

        Kmer::with_many_both(k, &seq, |x, y| {
            let c = CanonicalHash::new(k);
            let x_can = c.canonical(x);
            let y_can = c.canonical(y);
            assert_eq!(x_can, y_can);
        });
    }

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
        let m = (1u64 << (2 * k)) - 1;
        for _i in 0..20 {
            let u = rng.rnd() & m;
            let x = Kmer::from_u64(u);
            s += &x.render(k);
        }
        s
    }


    /*

    This one is expected to fail - the whole point of using
    hash based canonicalisation is to even out the distribution
    where you get a skewed distribution from lexicographic
    canonicalisation.

    #[test]
    fn test_lex_2() {
        let k = 15;

        let mut xs = Vec::new();
        let mut rng = MiniRng::new(19);
        let c = CanonicalLex::new(k);
        while xs.len() < 1024*1024 {
            let seq = make_some_sequence(&mut rng);
            Kmer::with_many(k, &seq, |x| {
                let y = c.canonical(x);
                xs.push(y);
            });
        }
        xs.sort();
        xs.dedup();
        let mut n = 0;
        let mut sx = 0;
        let mut sx2 = 0;
        for i in 1..xs.len() {
            let d = xs[i].0 - xs[i-1].0;
            n += 1;
            sx += d;
            sx2 += d*d;
        }
        let m = sx / n;
        let v = sx2 / n - m*m;
        let sd = (v as f64).sqrt();
        println!("mean = {}", m);
        println!("std.dev = {}", sd);
        assert_eq!(m, 1024);
        assert!((sd - 1024.0).abs() < 10.0);

    }
    */

    #[test]
    fn test_hash_2() {
        let k = 15;

        let mut xs = Vec::new();
        let mut rng = MiniRng::new(19);
        let c = CanonicalHash::new(k);
        while xs.len() < 1024*1024 {
            let seq = make_some_sequence(&mut rng);
            Kmer::with_many(k, &seq, |x| {
                let y = c.canonical(x);
                xs.push(y);
            });
        }
        xs.sort();
        xs.dedup();
        let mut n = 0;
        let mut sx = 0;
        let mut sx2 = 0;
        for i in 1..xs.len() {
            let d = xs[i].0 - xs[i-1].0;
            n += 1;
            sx += d;
            sx2 += d*d;
        }
        let m = sx / n;
        let v = sx2 / n - m*m;
        let sd = (v as f64).sqrt();
        assert_eq!(m, 1024);
        assert!((sd - 1024.0).abs() < 10.0);
    }
}