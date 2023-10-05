use crate::Kmer;

/// An adaptor that coverts an iterator over k-mers into a k-mer frequency iterator
pub fn unique<'a, Src>(itr: Src) -> Unique<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    Unique::new(itr)
}

/// An iterator producing k-mer frequencies.
pub struct Unique<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    src: Src,
}

impl<'a, Src> Unique<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    fn new(src: Src) -> Unique<'a, Src> {
        Unique { src }
    }
}

impl<'a, Src> Iterator for Unique<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    type Item = (Kmer, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.src.next().map(|x| (x.clone(), 1))
    }
}

/// Construct a k-mer frequency iterator from an iterator over frequencies.
pub fn frequency_vector_iter<'a, Src>(itr: Src) -> KmerFrequencyIterator<'a, Src>
where
    Src: Iterator<Item = &'a usize>,
{
    KmerFrequencyIterator::new(itr)
}

/// An iterator producing k-mer frequencies.
pub struct KmerFrequencyIterator<'a, Src>
where
    Src: Iterator<Item = &'a usize>,
{
    x: usize,
    src: Src,
}

impl<'a, Src> KmerFrequencyIterator<'a, Src>
where
    Src: Iterator<Item = &'a usize>,
{
    fn new(src: Src) -> KmerFrequencyIterator<'a, Src> {
        KmerFrequencyIterator { x: 0, src }
    }
}

impl<'a, Src> Iterator for KmerFrequencyIterator<'a, Src>
where
    Src: Iterator<Item = &'a usize>,
{
    type Item = (Kmer, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(freq) = self.src.next() {
            if *freq == 0 {
                self.x += 1;
                continue;
            }
            let x = Kmer::from_u64(self.x as u64);
            self.x += 1;
            return Some((x, *freq));
        }
        None
    }
}

/// Take a sorted iterator and yield k-mer frequencies.
pub fn counting_kmer_frequency_iterator<'a, Src>(src: Src) -> CountingIterator<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    CountingIterator::new(src)
}

/// Take a sorted iterator and yield k-mer frequencies.
pub struct CountingIterator<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    src: Src,
    src_next: Option<&'a Kmer>
}

impl<'a, Src> CountingIterator<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    fn new(mut src: Src) -> CountingIterator<'a, Src> {
        let src_next = src.next();
        CountingIterator { src, src_next }
    }
}

impl<'a, Src> Iterator for CountingIterator<'a, Src>
where
    Src: Iterator<Item = &'a Kmer>,
{
    type Item = (Kmer, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(x) = self.src_next {
            let x = x.clone();
            let mut f = 1;
            self.src_next = self.src.next();
            while let Some(y) = self.src_next {
                if *y != x {
                    break;
                }
                f += 1;
                self.src_next = self.src.next();
            }
            return Some((x, f));
        }
        None
    }
}

/// Merge two iterators
pub fn merge<Lhs, Rhs>(lhs: Lhs, rhs: Rhs) -> MergeIterator<Lhs, Rhs>
where
    Lhs: Iterator<Item = (Kmer, usize)>,
    Rhs: Iterator<Item = (Kmer, usize)>,
{
    MergeIterator::new(lhs, rhs)
}

/// Merge two k-mer frequency iterators.
pub struct MergeIterator<Lhs, Rhs>
where
    Lhs: Iterator<Item = (Kmer, usize)>,
    Rhs: Iterator<Item = (Kmer, usize)>,
{
    lhs: Lhs,
    lhs_next: Option<(Kmer, usize)>,
    rhs: Rhs,
    rhs_next: Option<(Kmer, usize)>,
}

impl<Lhs, Rhs> MergeIterator<Lhs, Rhs>
where
    Lhs: Iterator<Item = (Kmer, usize)>,
    Rhs: Iterator<Item = (Kmer, usize)>,
{
    fn new(mut lhs: Lhs, mut rhs: Rhs) -> MergeIterator<Lhs, Rhs> {
        let lhs_next = lhs.next();
        let rhs_next = rhs.next();
        MergeIterator {
            lhs,
            lhs_next,
            rhs,
            rhs_next,
        }
    }
}

impl<Lhs, Rhs> Iterator for MergeIterator<Lhs, Rhs>
where
    Lhs: Iterator<Item = (Kmer, usize)>,
    Rhs: Iterator<Item = (Kmer, usize)>,
{
    type Item = (Kmer, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some((x, x_f)), Some((y, y_f))) = (&self.lhs_next, &self.rhs_next) {
            if x < y {
                let res = (x.clone(), *x_f);
                self.lhs_next = self.lhs.next();
                return Some(res);
            }
            if y < x {
                let res = (y.clone(), *y_f);
                self.rhs_next = self.rhs.next();
                return Some(res);
            }
            let res = (x.clone(), *x_f + *y_f);
            self.lhs_next = self.lhs.next();
            self.rhs_next = self.rhs.next();
            return Some(res);
        }
        while let Some((x, x_f)) = &self.lhs_next {
            let res = (x.clone(), *x_f);
            self.lhs_next = self.lhs.next();
            return Some(res);
        }
        while let Some((y, y_f)) = &self.rhs_next {
            let res = (y.clone(), *y_f);
            self.rhs_next = self.rhs.next();
            return Some(res);
        }
        None
    }
}
