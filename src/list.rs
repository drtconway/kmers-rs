use std::io::Read;

use crate::{
    vbyte::{Decoder8, Encoder8},
    Kmer,
};

/// A compressed representation for a sorted list of k-mers.
pub struct CompressedKmerList {
    bytes: Vec<u8>,
    count: usize,
    last: Kmer,
}

impl CompressedKmerList {
    /// Create an empty CompressedKmerList
    pub fn new() -> CompressedKmerList {
        CompressedKmerList {
            bytes: Vec::new(),
            count: 0,
            last: Kmer(0),
        }
    }

    /// Return the number of k-mers in the list
    pub fn len(&self) -> usize {
        self.count
    }

    /// Add a k-mer to the compressed list.
    ///
    /// The k-mer `x` must be greater than or equal to the last k-mer on the list.
    pub fn push(&mut self, x: Kmer) {
        assert!(x >= self.last);
        let d = x.0 - self.last.0;
        Encoder8::encode(d, &mut self.bytes).unwrap();
        self.count += 1;
        self.last = x;
    }

    /// Create an iterator over the compressed list.
    pub fn iter(&self) -> CompressedKmerListIterator<'_, impl Iterator<Item = &'_ u8>> {
        CompressedKmerListIterator::new(self.bytes.iter())
    }

    /// Create a new k-mer list from an iterator.
    pub fn from_iter<Src>(src: Src) -> CompressedKmerList
    where
        Src: Iterator<Item = Kmer>,
    {
        let mut res = CompressedKmerList::new();
        for x in src {
            res.push(x);
        }
        res
    }
}

/// Iterate over a compressed k-mer list
pub struct CompressedKmerListIterator<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    bytes: IteratorReadAdaptor<'a, Src>,
    last: Kmer,
}

impl<'a, Src> CompressedKmerListIterator<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    /// Create an iterator over the compressed list.
    pub fn new(bytes: Src) -> CompressedKmerListIterator<'a, Src> {
        CompressedKmerListIterator {
            bytes: IteratorReadAdaptor::new(bytes),
            last: Kmer(0),
        }
    }
}

impl<'a, Src> Iterator for CompressedKmerListIterator<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    type Item = Kmer;
    fn next(&mut self) -> Option<Self::Item> {
        match Decoder8::decode(&mut self.bytes).unwrap() {
            None => None,
            Some(d) => {
                self.last.0 += d;
                Some(self.last.clone())
            }
        }
    }
}

/// A compressed representation for a sorted list of k-mer and count pairs.
pub struct CompressedKmerFrequencyList {
    bytes: Vec<u8>,
    count: usize,
    last: Kmer,
}

impl CompressedKmerFrequencyList {
    /// Create an empty CompressedKmerFrequencyList
    pub fn new() -> CompressedKmerFrequencyList {
        CompressedKmerFrequencyList {
            bytes: Vec::new(),
            count: 0,
            last: Kmer(0),
        }
    }

    /// Return the number of k-mer, frequency pairs in the list
    pub fn len(&self) -> usize {
        self.count
    }

    /// Add a k-mer and count to the compressed list.
    ///
    /// The k-mer `x` must be greater than or equal to the last k-mer on the list.
    pub fn push(&mut self, x_f: (Kmer, usize)) {
        let (x, f) = x_f;
        assert!(x >= self.last);
        let d = x.0 - self.last.0;
        Encoder8::encode(d, &mut self.bytes).unwrap();
        Encoder8::encode(f as u64, &mut self.bytes).unwrap();
        self.count += 1;
        self.last = x;
    }

    /// Create an iterator over the compressed list.
    pub fn iter(&self) -> CompressedKmerFrequencyListIterator<'_, impl Iterator<Item = &'_ u8>> {
        CompressedKmerFrequencyListIterator::new(self.bytes.iter())
    }

    /// Create a new k-mer list from an iterator.
    pub fn from_iter<Src>(src: Src) -> CompressedKmerFrequencyList
    where
        Src: Iterator<Item = (Kmer, usize)>,
    {
        let mut res = CompressedKmerFrequencyList::new();
        for x_f in src {
            res.push(x_f);
        }
        res
    }
}

/// Iterate over a compressed k-mer list
pub struct CompressedKmerFrequencyListIterator<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    bytes: IteratorReadAdaptor<'a, Src>,
    last: Kmer,
}

impl<'a, Src> CompressedKmerFrequencyListIterator<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    /// Create an iterator over the compressed list.
    pub fn new(bytes: Src) -> CompressedKmerFrequencyListIterator<'a, Src> {
        CompressedKmerFrequencyListIterator {
            bytes: IteratorReadAdaptor::new(bytes),
            last: Kmer(0),
        }
    }
}

impl<'a, Src> Iterator for CompressedKmerFrequencyListIterator<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    type Item = (Kmer, usize);
    fn next(&mut self) -> Option<Self::Item> {
        match Decoder8::decode(&mut self.bytes).unwrap() {
            None => None,
            Some(d) => {
                self.last.0 += d;
                let f = Decoder8::decode(&mut self.bytes).unwrap().unwrap();
                Some((self.last.clone(), f as usize))
            }
        }
    }
}

/// Convert an iterator over bytes to a [`Read`](std::io::Read).
pub struct IteratorReadAdaptor<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    itr: Src,
}

impl<'a, Src> IteratorReadAdaptor<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    /// Convert an iterator over bytes to a [`Read`](std::io::Read).
    pub fn new(src: Src) -> IteratorReadAdaptor<'a, Src> {
        IteratorReadAdaptor { itr: src }
    }
}

impl<'a, Src> Read for IteratorReadAdaptor<'a, Src>
where
    Src: Iterator<Item = &'a u8>,
{
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if buf.len() == 0 {
            return Ok(0);
        }

        let mut count = 0;
        while let Some(b) = self.itr.next() {
            buf[count] = *b;
            count += 1;
            if count == buf.len() {
                break;
            }
        }
        Ok(count)
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

        fn freq(&mut self) -> u64 {
            let mut f: f64 = 0.0;
            let mut u = self.rnd();
            while u > 0 {
                if u & 1 == 1 {
                    f += 1.0;
                }
                f /= 2.0;
                u >>= 1;
            }
            (1.0 / f).floor() as u64
        }
    }

    #[test]
    fn test_1() {
        let k = 7;
        let m = (1u64 << (2 * k)) - 1;
        let mut rng = MiniRng::new(19);
        let mut xs = Vec::new();
        for _i in 0..100 {
            let x = Kmer(rng.rnd() & m);
            xs.push(x);
        }
        xs.sort();
        let c = CompressedKmerList::from_iter(xs.iter().map(|x| x.clone()));
        assert_eq!(c.len(), xs.len());
        let ys = Vec::from_iter(c.iter());
        assert_eq!(ys, xs);
    }

    #[test]
    fn test_2() {
        let k = 7;
        let m = (1u64 << (2 * k)) - 1;
        let mut rng = MiniRng::new(19);
        let mut xs = Vec::new();
        for _i in 0..100 {
            let x = Kmer(rng.rnd() & m);
            let f = rng.freq();
            xs.push((x, f as usize));
        }
        xs.sort();
        let c = CompressedKmerFrequencyList::from_iter(xs.iter().map(|(x,f)| (x.clone(), *f)));
        assert_eq!(c.len(), xs.len());
        let ys = Vec::from_iter(c.iter());
        assert_eq!(ys, xs);
    }}
