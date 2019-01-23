pub mod counts;
pub mod features;

use std::collections::BTreeMap;

use self::{
    counts::Counts,
    features::{merge_intervals, Feature, Features},
};

pub fn calculate_fpkms(
    counts: &Counts,
    features: &Features,
    scaling_factor: f64,
) -> BTreeMap<String, f64> {
    counts
        .iter()
        .map(|(name, &count)| {
            let intervals = &features[name];
            let len = sum_nonoverlapping_interval_lengths(intervals);
            let fpkm = calculate_fpkm(count, len, scaling_factor);
            (name.clone(), fpkm)
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

    use super::*;

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
