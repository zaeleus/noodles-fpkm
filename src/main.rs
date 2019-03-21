use std::{
    fs::File,
    io::{self, Write},
};

use clap::{crate_name, crate_version, App, Arg};
use log::LevelFilter;
use noodles_fpkm::{
    calculate_fpkms, calculate_tpms, counts::read_counts, features::read_features, Expressions,
};

fn write_expressions<W>(mut writer: W, expressions: &Expressions) -> io::Result<()>
where
    W: Write,
{
    for (id, value) in expressions {
        writeln!(writer, "{}\t{}", id, value)?;
    }

    Ok(())
}

fn main() {
    let matches = App::new(crate_name!())
        .version(crate_version!())
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Use verbose logging"),
        )
        .arg(
            Arg::with_name("feature-type")
                .short("t")
                .long("type")
                .value_name("str")
                .help("Feature type to count")
                .default_value("exon"),
        )
        .arg(
            Arg::with_name("feature-id")
                .short("i")
                .long("id")
                .value_name("str")
                .help("Feature attribute to use as the feature identity")
                .default_value("gene_id"),
        )
        .arg(
            Arg::with_name("method")
                .short("m")
                .long("method")
                .value_name("str")
                .help("Normalization method")
                .default_value("fpkm")
                .possible_values(&["fpkm", "tpm"]),
        )
        .arg(
            Arg::with_name("annotations")
                .short("a")
                .long("annotations")
                .value_name("file")
                .help("Input annotations file (GTF/GFFv2)")
                .required(true),
        )
        .arg(
            Arg::with_name("counts")
                .help("Input feature counts")
                .required(true)
                .index(1),
        )
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
    let unit = matches.value_of("unit").unwrap();

    let features = read_features(annotations_src, feature_type, feature_id).unwrap();

    let file = File::open(&counts_src).unwrap();
    let counts = read_counts(file).unwrap();

    let fpkms = match unit {
        "fpkm" => calculate_fpkms(&counts, &features).unwrap(),
        "tpm" => calculate_tpms(&counts, &features).unwrap(),
        _ => unreachable!(),
    };

    let stdout = io::stdout();
    let handle = stdout.lock();
    write_expressions(handle, &fpkms).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_expressions() {
        let fpkms = [
            (String::from("AAAS"), 5825.440538780093),
            (String::from("AC009952.3"), 10.494073576888187),
            (String::from("RPL37AP1"), 3220170.8708099453),
            (String::from("ZNF700"), 0.0),
        ]
        .iter()
        .cloned()
        .collect();

        let mut buf = Vec::new();
        write_expressions(&mut buf, &fpkms).unwrap();

        let actual = String::from_utf8(buf).unwrap();
        let expected = "\
AAAS\t5825.440538780093
AC009952.3\t10.494073576888187
RPL37AP1\t3220170.8708099453
ZNF700\t0
";

        assert_eq!(actual, expected);
    }
}
