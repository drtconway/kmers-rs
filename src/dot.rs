use crate::Kmer;

/// Compute the dot product between two k-mer frequency spectra.
pub fn dot<Lhs, Rhs>(mut lhs: Lhs, mut rhs: Rhs) -> f64
where
    Lhs: Iterator<Item = (Kmer, usize)>,
    Rhs: Iterator<Item = (Kmer, usize)>,
{
    let mut d = 0.0;

    let mut lhs_next = lhs.next();
    let mut rhs_next = rhs.next();
    while let (Some(x), Some(y)) = (&lhs_next, &rhs_next) {
        if x.0 < y.0 {
            lhs_next = lhs.next();
        } else if x.0 > y.0 {
            rhs_next = rhs.next();
        } else {
            d += (x.1 * y.1) as f64;
            lhs_next = lhs.next();
            rhs_next = rhs.next();
        }
    }

    d
}

#[cfg(test)]
mod test {
    use crate::frequency::{unique, frequency_vector_iter};

    use super::*;

    #[test]
    fn test_1() {
        let k = 11;

        let lhs_seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
        let mut lhs_kmers = Kmer::make_many(k, lhs_seq);
        lhs_kmers.sort();
        lhs_kmers.dedup();

        let rhs_seq = "CCACGACACATCACAGTCAAATTTCTAAAAATTAAGGACCACACAAACGCCTTAAAAACAGAGAAAGAAAAATGACATATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGCAACCACAGAGTTCAGAATGAAG";
        let mut rhs_kmers = Kmer::make_many(k, rhs_seq);
        rhs_kmers.sort();
        rhs_kmers.dedup();

        let lhs_kmer_frequencies = unique(lhs_kmers.iter().map(|x| x.clone()));
        let rhs_kmer_frequencies = unique(rhs_kmers.iter().map(|x| x.clone()));
        let d = dot(lhs_kmer_frequencies, rhs_kmer_frequencies);
        assert_eq!(d, 96.0);
    }

    #[test]
    fn test_2() {
        let k = 5;

        let lhs_seq = "CCACGACACATCATAGTCAAATTTCTAAAAATTAAGGACCACACAAACACCTTAAAAACAGAGAAAGAAAAATGACACATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGAAACCACAGAGTTCAGAATGAAG";
        let lhs_frequencies = Kmer::frequency_vector(k, &lhs_seq);

        let rhs_seq = "CCACGACACATCACAGTCAAATTTCTAAAAATTAAGGACCACACAAACGCCTTAAAAACAGAGAAAGAAAAATGACATATTACCTATACTGCAAAAACAATTCAAATGACAAAGTATTATTCATTGGCAACCACAGAGTTCAGAATGAAG";
        let rhs_frequencies = Kmer::frequency_vector(k, &rhs_seq);

        let lhs_kmer_frequencies = frequency_vector_iter(lhs_frequencies.iter());
        let rhs_kmer_frequencies = frequency_vector_iter(rhs_frequencies.iter());
        let d = dot(lhs_kmer_frequencies, rhs_kmer_frequencies);
        assert_eq!(d, 188.0);
    }
}
