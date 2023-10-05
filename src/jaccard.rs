use crate::Kmer;

/// Compute the Jaccard cooeficient from two sets of k-mers.
pub fn jaccard<'a, Lhs, Rhs>(mut lhs: Lhs, mut rhs: Rhs) -> f64
where
    Lhs: Iterator<Item = &'a Kmer>,
    Rhs: Iterator<Item = &'a Kmer>,
{
    let mut i = 0;
    let mut u = 0;

    let mut lhs_next = lhs.next();
    let mut rhs_next = rhs.next();
    while let (Some(x), Some(y)) = (&lhs_next, &rhs_next) {
        if x < y {
            u += 1;
            lhs_next = lhs.next();
        } else if x > y {
            u += 1;
            rhs_next = rhs.next();
        } else {
            i += 1;
            u += 1;
            lhs_next = lhs.next();
            rhs_next = rhs.next();
        }
    }
    while let Some(_x) = lhs_next {
        u += 1;
        lhs_next = lhs.next();
    }
    while let Some(_y) = rhs_next {
        u += 1;
        rhs_next = rhs.next();
    }
    (i as f64) / (u as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::HashSet;

    #[test]
    fn test_1() {
        let k = 11;

        let lhs_seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
        let mut lhs_kmers = Kmer::make_many(k, lhs_seq);
        lhs_kmers.sort();
        lhs_kmers.dedup();
        let lhs_set: HashSet<Kmer> = HashSet::from_iter(lhs_kmers.iter().map(|x| x.clone()));

        let rhs_seq = "CCACGACACATCACAGTCAAATTTCTAAAAATTAAGGACCACACAAACGCCTTAAAAACAGAGAAAGAAAAATGACATATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGCAACCACAGAGTTCAGAATGAAG";
        let mut rhs_kmers = Kmer::make_many(k, rhs_seq);
        rhs_kmers.sort();
        rhs_kmers.dedup();
        let rhs_set: HashSet<Kmer> = HashSet::from_iter(rhs_kmers.iter().map(|x| x.clone()));

        let i = rhs_set.intersection(&lhs_set).count();
        let u = rhs_set.union(&lhs_set).count();
        let j = (i as f64) / (u as f64);

        let d = jaccard(lhs_kmers.iter(), rhs_kmers.iter());
        assert_eq!(d, j);
    }
}
