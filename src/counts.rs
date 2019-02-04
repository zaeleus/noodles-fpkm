use std::{
    collections::hash_map::{Entry, HashMap},
    io::{self, Read},
};

use csv::StringRecord;

const NAME_INDEX: usize = 0;
const COUNT_INDEX: usize = 1;

static HTSEQ_COUNT_META_PREFIX: &str = "__";

pub type Counts = HashMap<String, u64>;

/// Reads TSV-formatted data and returns a map of feature ID-count pairs.
///
/// The input is TSV-formatted with two columns: a feature identifier (string)
/// and a count (integer).
///
/// Reading stops at EOF or the first identifier that starts with "__". This
/// prefix is considered to be a special counter or extra metadata, as defined
/// by [htseq-count] > 0.5.4.
///
/// [htseq-count]: https://htseq.readthedocs.io/en/release_0.11.1/count.html#usage
///
/// # Example
///
/// ```
/// use noodles_fpkm::counts::read_counts;
///
/// let data = "\
/// AAAS\t645
/// AC009952.3\t1
/// RPL37AP1\t5714
/// __no_feature\t136550
/// ";
///
/// let counts = read_counts(data.as_bytes()).unwrap();
///
/// assert_eq!(counts.len(), 3);
/// assert_eq!(counts["AAAS"], 645);
/// assert_eq!(counts["AC009952.3"], 1);
/// assert_eq!(counts["RPL37AP1"], 5714);
/// ```
pub fn read_counts<R>(reader: R) -> io::Result<Counts> where R: Read {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(reader);

    let mut counts = Counts::new();

    for result in rdr.records() {
        let record = result?;

        let name = parse_name(&record)?;

        if name.starts_with(HTSEQ_COUNT_META_PREFIX) {
            break;
        }

        let count = parse_count(&record)?;

        insert_count(&mut counts, name, count)?;
    }

    Ok(counts)
}

fn parse_name(record: &StringRecord) -> io::Result<&str> {
    let cell = record.get(NAME_INDEX);

    cell.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid name: {:?}", cell),
        )
    })
}

fn parse_count(record: &StringRecord) -> io::Result<u64> {
    let cell = record.get(COUNT_INDEX);

    cell.and_then(|s| s.parse().ok()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid count: {:?}", cell),
        )
    })
}

fn insert_count<'a>(counts: &'a mut Counts, name: &str, count: u64) -> io::Result<&'a mut u64> {
    match counts.entry(name.to_string()) {
        Entry::Vacant(e) => Ok(e.insert(count)),
        Entry::Occupied(_) => {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("duplicate identifier '{}'", name),
            ))
        },
    }
}

/// Sums the counts from a `Count` map.
///
/// # Example
///
/// ```
/// use noodles_fpkm::counts::sum_counts;
///
/// let counts = [
///     (String::from("AAAS"), 645),
///     (String::from("AC009952.3"), 1),
///     (String::from("RPL37AP1"), 5714),
/// ].iter().cloned().collect();
///
/// assert_eq!(sum_counts(&counts), 6360);
/// ```
pub fn sum_counts(counts: &Counts) -> u64 {
    counts.values().sum()
}

#[cfg(test)]
mod tests {
    use csv::StringRecord;

    use super::*;

    #[test]
    fn test_read_counts_with_duplicate_identifiers() {
        let data = "\
AAAS\t645
AC009952.3\t1
AC009952.3\t0
RPL37AP1\t5714
__no_feature\t136550
";

        assert!(read_counts(data.as_bytes()).is_err());
    }

    #[test]
    fn test_parse_name() {
        let record = StringRecord::from(vec!["AAAS", "645"]);
        assert_eq!(parse_name(&record).unwrap(), "AAAS");
    }

    #[test]
    fn test_parse_count() {
        let record = StringRecord::from(vec!["AAAS", "645"]);
        assert_eq!(parse_count(&record).unwrap(), 645);

        let record = StringRecord::from(vec!["AAAS", ""]);
        assert!(parse_count(&record).is_err());

        let record = StringRecord::from(vec!["AAAS", "x"]);
        assert!(parse_count(&record).is_err());
    }

    #[test]
    fn test_insert_count() {
        let mut counts = Counts::new();

        assert!(insert_count(&mut counts, "AAAS", 645).is_ok());
        assert_eq!(counts["AAAS"], 645);

        assert!(insert_count(&mut counts, "AAAS", 0).is_err());
    }
}
