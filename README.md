**Deprecation notice**: noodles-fpkm was merged into [noodles-squab]. Use
`noodles-squab normalize` instead.

[noodles-squab]: https://github.com/zaeleus/noodles-squab

---

# noodles fpkm

**noodles fpkm** calculates FPKM (**f**ragments **p**er **k**ilobase per
**m**illion mapped reads) values from feature counts.

## Install

Install [Rust] and use `cargo` to install `noodles-fpkm`.

```
$ cargo install --git https://github.com/zaeleus/noodles-fpkm.git
```

[Rust]: https://www.rust-lang.org/tools/install

## Usage

```
noodles-fpkm 0.1.0

USAGE:
    noodles-fpkm [FLAGS] [OPTIONS] <counts> --annotations <file>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    Use verbose logging

OPTIONS:
    -a, --annotations <file>    Input annotations file (GTF/GFFv2)
    -i, --id <str>              Feature attribute to use as the feature identity [default: gene_id]
    -t, --type <str>            Feature type to count [default: exon]

ARGS:
    <counts>    Input feature counts
```

Output is printed to `stdout` as tab-separated values with no header. There are
two columns: the feature identifier (string) and its FPKM value (float), e.g.,

```
AAAS    3.568542663590036
AC009952.3      0.006428449320644941
AC009952.4      0
RPL37AP1        27.617908579562194
```

The table is sorted lexographically by the feature identifier.

## Example

Use [noodles-count-features] or [htseq-count] to create a table of feature
counts.

```
$ noodles-count-features \
  --annotations annotations.gtf \
  --type gene \
  --id gene_name \
  --output counts.txt \
  sample.bam
```

Use the same annotations (and options) with the resulting counts to calculate
FPKM values.

```
$ noodles-fpkm \
  --annotations annotations.gtf \
  --type gene \
  --id gene_name \
  counts.txt
```

[noodles-count-features]: https://github.com/zaeleus/noodles-count-features
[htseq-count]: https://htseq.readthedocs.io/en/release_0.11.1/count.html

## References

  * Wagner, G.P., Kin, K. & Lynch, V.J. Theory Biosci. (2012) 131: 281.
    https://doi.org/10.1007/s12064-012-0162-3
