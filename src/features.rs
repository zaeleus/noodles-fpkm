use std::{collections::HashMap, io, path::Path};

use noodles::formats::gff;

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

    /// Returns the size of the feature's interval.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_fpkm::features::Feature;
    ///
    /// assert_eq!(Feature::new(2, 5).len(), 4);
    /// assert_eq!(Feature::new(3, 4).len(), 2);
    /// assert_eq!(Feature::new(5, 7).len(), 3);
    /// assert_eq!(Feature::new(9, 12).len(), 4);
    /// assert_eq!(Feature::new(10, 15).len(), 6);
    /// assert_eq!(Feature::new(19, 21).len(), 3);
    /// ```
    pub fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}

/// Merges a list of overlapping intervals into a list of non-overlapping intervals.
///
/// The intervals are assumed to be inclusive.
///
/// # Panics
///
/// Panics when the input list of intervals is empty.
///
/// # Example
///
/// Given the list of intervals [2, 5], [3, 4], [5, 7], [9, 12], [10, 15], and
/// [16, 21], the following overlap and are combined into single intervals.
///
///
///   * [2, 5], [3, 4], [5, 7] => [2, 7]
///   * [9, 12], [10, 15] => [9, 15]
///   * [16, 21] => [16, 21]
///
/// ```
/// use noodles_fpkm::features::{merge_intervals, Feature};
///
/// let features = [
///     Feature::new(2, 5), Feature::new(3, 4), Feature::new(5, 7),
///     Feature::new(9, 12), Feature::new(10, 15), Feature::new(16, 21),
/// ];
///
/// let actual = merge_intervals(&features);
/// let expected = [Feature::new(2, 7), Feature::new(9, 15), Feature::new(16, 21)];
/// assert_eq!(actual, expected);
/// ```
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


/// Builds a map of feature ID-feature vector pairs from a GTF/GFFv2.
///
/// The [GTF/GFFv2] is filtered by `feature_type` (column 3), using
/// `feature_id` as the key for the map from the feature attributes
/// (column 9).
///
/// [GTF/GFFv2]: https://useast.ensembl.org/info/website/upload/gff.html
///
/// # Example
///
/// ```
/// use noodles_fpkm::features::{read_features, Feature};
///
/// let features = read_features(
///     "test/fixtures/annotations.gtf",
///     "exon",
///     "gene_name",
/// ).unwrap();
///
/// assert_eq!(features.len(), 2);
///
/// assert_eq!(
///     &features["DDX11L1"],
///     &[Feature::new(11869, 12227), Feature::new(12613, 12721)],
/// );
///
/// assert_eq!(&features["NECAP2"], &[Feature::new(16440672, 16440853)]);
/// ```
pub fn read_features<P>(src: P, feature_type: &str, feature_id: &str) -> io::Result<Features>
where
    P: AsRef<Path>,
{
    let mut reader = gff::open(src)?;
    let mut features: Features = HashMap::new();

    for result in reader.records() {
        let row = result?;
        let record = gff::Record::new(row);

        let ty = record.feature().map_err(invalid_data)?;

        if ty != feature_type {
            continue;
        }

        let start = record.start().map_err(invalid_data)?;
        let end = record.end().map_err(invalid_data)?;

        let attributes = record.attributes().map_err(invalid_data)?;
        let id = attributes.get(feature_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing attribute '{}'", feature_id),
            )
        })?;

        let list = features.entry(id.to_string()).or_default();
        let feature = Feature::new(start, end);
        list.push(feature);
    }

    Ok(features)
}

fn invalid_data(e: gff::record::Error) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e)
}
