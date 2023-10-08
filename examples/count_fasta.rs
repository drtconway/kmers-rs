extern crate noodles;

use std::collections::BTreeMap;

use kmers::BigSparseAccumulator;

fn main() -> std::io::Result<()> {
    let args = Vec::from_iter(std::env::args());
    if args.len() != 2 {
        println!("to run this example, download a FASTA file such as the bacterial genome at");
        println!("    https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP113114.1?report=fasta");
        println!("then invoke the code in the following manner:");
        println!("    cargo run --example count_fasta -- NZ_CP113114.1.fasta");
    }

    let k = 15;
    let buffer_size = 1024*1024;

    let mut reader = noodles::fasta::reader::Builder.build_from_path(&args[1])?;
    let mut acc = BigSparseAccumulator::new(k, buffer_size);
    for res in reader.records() {
        let rec = res?;
        acc.add_both(rec.sequence());
    }
    let mut hist = BTreeMap::new();
    for (_x,f) in acc.kmer_frequencies() {
        *hist.entry(f).or_insert(0) += 1;
    }
    for (f, c) in hist.iter() {
        println!("{}\t{}", f, c);
    }
    Ok(())
}