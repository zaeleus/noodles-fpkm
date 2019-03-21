pub mod counts;
pub mod features;

use std::collections::{BTreeMap, HashMap};

use self::{
    counts::{sum_counts, Counts},
    features::{merge_intervals, Feature, Features},
};

#[derive(Debug)]
pub enum Error {
    MissingFeature(String),
}

pub type Expressions = BTreeMap<String, f64>;

pub fn calculate_fpkms(counts: &Counts, features: &Features) -> Result<Expressions, Error> {
    let counts_sum = sum_counts(counts);

    counts
        .iter()
        .map(|(name, &count)| {
            features
                .get(name)
                .map(|intervals| {
                    let len = sum_nonoverlapping_interval_lengths(intervals);
                    let fpkm = calculate_fpkm(count, len, counts_sum);
                    (name.clone(), fpkm)
                })
                .ok_or_else(|| Error::MissingFeature(name.clone()))
        })
        .collect()
}

fn sum_nonoverlapping_interval_lengths(intervals: &[Feature]) -> u64 {
    merge_intervals(intervals).iter().map(|i| i.len()).sum()
}

// See <https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/>.
fn calculate_fpkm(count: u64, len: u64, counts_sum: u64) -> f64 {
    (count as f64 * 1e9) / (len as f64 * counts_sum as f64)
}

pub fn calculate_tpms(counts: &Counts, features: &Features) -> Result<Expressions, Error> {
    let cpbs: HashMap<String, f64> = counts
        .iter()
        .map(|(name, &count)| {
            features
                .get(name)
                .map(|intervals| {
                    let len = sum_nonoverlapping_interval_lengths(intervals);
                    let cpb = count as f64 / len as f64;
                    (name.clone(), cpb)
                })
                .ok_or_else(|| Error::MissingFeature(name.clone()))
        })
        .collect::<Result<_, _>>()?;

    let cpbs_sum = cpbs.values().sum();

    let tpms = cpbs
        .iter()
        .map(|(name, &cpb)| (name.clone(), calculate_tpm(cpb, cpbs_sum)))
        .collect();

    Ok(tpms)
}

fn calculate_tpm(cpb: f64, cpbs_sum: f64) -> f64 {
    cpb * 1e6 / cpbs_sum
}

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use crate::{
        counts::Counts,
        features::{Feature, Features},
    };

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
            (
                String::from("AC009952.3"),
                vec![Feature::new(9189629, 9204611)],
            ),
            (
                String::from("RPL37AP1"),
                vec![Feature::new(44466564, 44466842)],
            ),
        ];

        features.into_iter().cloned().collect()
    }

    #[test]
    fn test_calculate_fpkms() {
        let counts = build_counts();
        let features = build_features();

        let fpkms = calculate_fpkms(&counts, &features).unwrap();

        assert_eq!(fpkms.len(), 3);

        let a = fpkms["AAAS"];
        let b = 5825.440538780093;
        assert!((a - b).abs() < EPSILON);

        let a = fpkms["AC009952.3"];
        let b = 10.494073576888189;
        assert!((a - b).abs() < EPSILON);

        let a = fpkms["RPL37AP1"];
        let b = 3220170.8708099457;
        assert!((a - b).abs() < EPSILON);
    }

    #[test]
    fn test_calculate_fpkms_is_ordered_by_feature_id() {
        let counts = build_counts();
        let features = build_features();

        let fpkms = calculate_fpkms(&counts, &features).unwrap();
        let mut ids = fpkms.keys();

        assert_eq!(ids.next().unwrap(), "AAAS");
        assert_eq!(ids.next().unwrap(), "AC009952.3");
        assert_eq!(ids.next().unwrap(), "RPL37AP1");
        assert!(ids.next().is_none());
    }

    #[test]
    fn test_calculate_fpkms_with_missing_feature() {
        let counts = build_counts();

        let mut features = build_features();
        features.remove("AC009952.3");

        assert!(calculate_fpkms(&counts, &features).is_err());
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
        let a = calculate_fpkm(2, 10, 212);
        let b = 943396.2264150943;
        assert!((a - b).abs() < EPSILON);

        let a = calculate_fpkm(5, 138756, 600081);
        let b = 0.06004935631747696;
        assert!((a - b).abs() < EPSILON);
    }

    #[test]
    fn test_calculate_tpm() {
        let a = calculate_tpm(2.0, 10.0);
        let b = 200000.0;
        assert!((a - b).abs() < EPSILON);

        let a = dbg!(calculate_tpm(0.0010, 26.65));
        let b = 37.5234521575985;
        assert!((a - b).abs() < EPSILON);
    }
}
