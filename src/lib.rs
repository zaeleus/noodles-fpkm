pub mod counts;
pub mod features;

use std::collections::BTreeMap;

use self::{
    counts::Counts,
    features::{merge_intervals, Feature, Features},
};

#[derive(Debug)]
pub enum Error {
    MissingFeature(String),
}

pub fn calculate_fpkms(
    counts: &Counts,
    features: &Features,
    scaling_factor: f64,
) -> Result<BTreeMap<String, f64>, Error> {
    counts
        .iter()
        .map(|(name, &count)| {
            features
                .get(name)
                .map(|intervals| {
                    let len = sum_nonoverlapping_interval_lengths(intervals);
                    let fpkm = calculate_fpkm(count, len, scaling_factor);
                    (name.clone(), fpkm)
                })
                .ok_or_else(|| Error::MissingFeature(name.clone()))
        })
        .collect()
}

fn sum_nonoverlapping_interval_lengths(intervals: &[Feature]) -> u64 {
    merge_intervals(intervals).iter().map(|i| i.len()).sum()
}

// See <https://statquest.org/2015/07/09/rpkm-fpkm-and-tpm-clearly-explained/>.
fn calculate_fpkm(count: u64, len: u64, scaling_factor: f64) -> f64 {
    let rpm = (count as f64) / scaling_factor;
    let len_kb = (len as f64) / 1000.0;
    rpm / len_kb
}

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use crate::{counts::Counts, features::{Feature, Features}};

    use super::*;

    fn build_counts() -> Counts {
        let counts = [
            (String::from("AAAS"), 645),
            (String::from("AC009952.3"), 1),
            (String::from("RPL37AP1"), 5714),
        ];

        counts.iter().cloned().collect()
    }

    fn build_features() -> Features {
        let features = [
            (String::from("AAAS"), vec![Feature::new(53307456, 53324864)]),
            (String::from("AC009952.3"), vec![Feature::new(9189629, 9204611)]),
            (String::from("RPL37AP1"), vec![Feature::new(44466564, 44466842)]),
        ];

        features.into_iter().cloned().collect()
    }

    #[test]
    fn test_calculate_fpkms() {
        let counts = build_counts();
        let features = build_features();
        let scaling_factor = 6360.0 / 1_000_000.0;

        let fpkms = calculate_fpkms(&counts, &features, scaling_factor).unwrap();

        assert_eq!(fpkms.len(), 3);

        let a = 5825.440538780093;
        let b = fpkms["AAAS"];
        assert!((a - b).abs() < EPSILON);

        let a = 10.494073576888187;
        let b = fpkms["AC009952.3"];
        assert!((a - b).abs() < EPSILON);

        let a = 3220170.8708099453;
        let b = fpkms["RPL37AP1"];
        assert!((a - b).abs() < EPSILON);
    }

    #[test]
    fn test_calculate_fpkms_with_missing_feature() {
        let counts = build_counts();

        let mut features = build_features();
        features.remove("AC009952.3");

        let scaling_factor = 6360.0 / 1_000_000.0;

        assert!(calculate_fpkms(&counts, &features, scaling_factor).is_err());
    }

    #[test]
    fn test_sum_nonoverlapping_interval_lengths() {
        let features = [
            Feature::new(2, 5),
            Feature::new(3, 4),
            Feature::new(5, 7),
            Feature::new(9, 12),
            Feature::new(10, 15),
            Feature::new(16, 21),
        ];

        let len = sum_nonoverlapping_interval_lengths(&features);
        assert_eq!(len, 19);
    }

    #[test]
    fn test_calculate_fpkm() {
        let a = calculate_fpkm(2, 10, 2.0);
        let b = 100.0;
        assert!((a - b).abs() < EPSILON);

        let a = calculate_fpkm(2, 4364, 10.382334);
        let b = 0.04414182225995562;
        assert!((a - b).abs() < EPSILON);
    }
}
