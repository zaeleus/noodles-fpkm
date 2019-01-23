use std::{collections::HashMap, io, path::Path};

use csv::StringRecord;
use noodles::formats::gff;

const GFF_GENE_INDEX: usize = 2;
const GFF_START_INDEX: usize = 3;
const GFF_END_INDEX: usize = 4;
const GFF_ATTRS_INDEX: usize = 8;

pub type Features = HashMap<String, Vec<Feature>>;

// 1-based, inclusive
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Feature {
    pub start: u64,
    pub end: u64,
}

impl Feature {
    pub fn new(start: u64, end: u64) -> Feature {
        Feature { start, end }
    }

    pub fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}

#[cfg(test)]
mod feature_tests {
    use super::Feature;

    #[test]
    fn test_len() {
        assert_eq!(Feature::new(2, 5).len(), 4);
        assert_eq!(Feature::new(3, 4).len(), 2);
        assert_eq!(Feature::new(5, 7).len(), 3);
        assert_eq!(Feature::new(9, 12).len(), 4);
        assert_eq!(Feature::new(10, 15).len(), 6);
        assert_eq!(Feature::new(19, 21).len(), 3);
    }
}

pub fn merge_intervals(intervals: &[Feature]) -> Vec<Feature> {
    assert!(!intervals.is_empty());

    let mut intervals = intervals.to_vec();
    intervals.sort_unstable_by_key(|i| i.start);

    let mut merged_intervals = Vec::with_capacity(intervals.len());
    merged_intervals.push(intervals[0].clone());

    for b in intervals {
        let a = merged_intervals.last_mut().expect("list cannot be empty");

        if b.start > a.end {
            merged_intervals.push(b);
            continue;
        }

        if a.end < b.end {
            a.end = b.end;
        }
    }

    merged_intervals
}

pub fn read_features<P>(src: P, feature_type: &str, feature_id: &str) -> io::Result<Features>
where
    P: AsRef<Path>,
{
    let mut reader = gff::open(src)?;
    let mut features: Features = HashMap::new();

    for result in reader.records() {
        let record = result?;

        let ty = &record[GFF_GENE_INDEX];

        if ty != feature_type {
            continue;
        }

        let start = parse_start(&record)?;
        let end = parse_end(&record)?;
        let id = parse_attrs_and_get(&record, feature_id)?;

        let list = features.entry(id.to_string()).or_default();
        let feature = Feature::new(start, end);
        list.push(feature);
    }

    Ok(features)
}

fn parse_start(record: &StringRecord) -> io::Result<u64> {
    let cell = record.get(GFF_START_INDEX);

    cell.and_then(|s| s.parse().ok()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid start: {:?}", cell),
        )
    })
}

fn parse_end(record: &StringRecord) -> io::Result<u64> {
    let cell = record.get(GFF_END_INDEX);

    cell.and_then(|s| s.parse().ok()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid end: {:?}", cell),
        )
    })
}

fn parse_attrs_and_get<'a>(record: &'a StringRecord, key: &str) -> io::Result<&'a str> {
    let attrs = record
        .get(GFF_ATTRS_INDEX)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid attrs"))?;

    for attr in attrs.split(';').map(str::trim_left) {
        if attr.starts_with(key) {
            return attr.split('"').nth(1).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("could not parse attribute '{}'", attr),
                )
            });
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("missing attribute '{}'", key),
    ))
}

#[cfg(test)]
mod tests {
    use csv::StringRecord;

    use super::*;

    #[test]
    fn test_merge_intervals() {
        let features = [
            Feature::new(2, 5),
            Feature::new(3, 4),
            Feature::new(5, 7),
            Feature::new(9, 12),
            Feature::new(10, 15),
            Feature::new(16, 21),
        ];

        let actual = merge_intervals(&features);

        let expected = [
            Feature::new(2, 7),
            Feature::new(9, 15),
            Feature::new(16, 21),
        ];

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_read_features() {
        let features = read_features(
            "test/fixtures/annotations.gtf",
            "exon",
            "gene_name",
        ).unwrap();

        assert_eq!(features.len(), 2);

        assert_eq!(
            &features["DDX11L1"],
            &[Feature::new(11869, 12227), Feature::new(12613, 12721)],
        );

        assert_eq!(&features["NECAP2"], &[Feature::new(16440672, 16440853)]);
    }

    fn build_record() -> StringRecord {
        StringRecord::from(vec![
           "chr1",
           "HAVANA",
           "gene",
           "11869",
           "14409",
           ".",
           "+",
           ".",
           r#"gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"#
        ])
    }

    #[test]
    fn test_parse_start() {
        let record = build_record();
        assert_eq!(parse_start(&record).unwrap(), 11869);
    }

    #[test]
    fn test_parse_end() {
        let record = build_record();
        assert_eq!(parse_end(&record).unwrap(), 14409);
    }

    #[test]
    fn test_parse_attrs_and_get() {
        let record = build_record();

        assert_eq!(parse_attrs_and_get(&record, "gene_id").unwrap(), "ENSG00000223972.5");
        assert_eq!(parse_attrs_and_get(&record, "gene_type").unwrap(), "transcribed_unprocessed_pseudogene");
        assert_eq!(parse_attrs_and_get(&record, "gene_name").unwrap(), "DDX11L1");
        assert_eq!(parse_attrs_and_get(&record, "havana_gene").unwrap(), "OTTHUMG00000000961.2");

        assert!(parse_attrs_and_get(&record, "level").is_err());
        assert!(parse_attrs_and_get(&record, "dne").is_err());
    }
}
