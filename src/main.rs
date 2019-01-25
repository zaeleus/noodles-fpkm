use std::fs::File;

use clap::{crate_name, crate_version, App, Arg};
use log::LevelFilter;
use noodles_fpkm::{
    calculate_fpkms,
    counts::{read_counts, sum_counts},
    features::read_features,
};

fn main() {
    let matches = App::new(crate_name!())
        .version(crate_version!())
        .arg(Arg::with_name("verbose")
             .short("v")
             .long("verbose")
             .help("Use verbose logging"))
        .arg(Arg::with_name("feature-type")
             .short("t")
             .long("type")
             .value_name("str")
             .help("Feature type to count")
             .default_value("exon"))
        .arg(Arg::with_name("feature-id")
             .short("i")
             .long("id")
             .value_name("str")
             .help("Feature attribute to use as the feature identity")
             .default_value("gene_id"))
        .arg(Arg::with_name("annotations")
             .short("a")
             .long("annotations")
             .value_name("file")
             .help("Input annotations file (GTF/GFFv2)")
             .required(true))
        .arg(Arg::with_name("counts")
             .help("Input feature counts")
             .required(true)
             .index(1))
        .get_matches();

    if matches.is_present("verbose") {
        env_logger::Builder::from_default_env()
            .filter(Some("noodles_fpkm"), LevelFilter::Info)
            .init();
    } else {
        env_logger::init();
    }

    let counts_src = matches.value_of("counts").unwrap();
    let annotations_src = matches.value_of("annotations").unwrap();
    let feature_type = matches.value_of("feature-type").unwrap();
    let feature_id = matches.value_of("feature-id").unwrap();

    let features = read_features(annotations_src, feature_type, feature_id).unwrap();

    let file = File::open(&counts_src).unwrap();
    let counts = read_counts(file).unwrap();

    let counts_sum = sum_counts(&counts) as f64;
    let scaling_factor = counts_sum / 1_000_000.0;

    let fpkms = calculate_fpkms(&counts, &features, scaling_factor).unwrap();

    for (name, fpkm) in fpkms {
        println!("{}\t{}", name, fpkm);
    }
}
