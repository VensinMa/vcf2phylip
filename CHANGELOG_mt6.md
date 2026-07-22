# 2.9-mt6 changes

- Detect plain VCF, ordinary gzip, and BGZF from file bytes rather than filename suffix.
- Validate the BGZF `BC` extra subfield and first compressed block.
- Read plain VCF as binary chunks and transfer contiguous byte buffers to parser workers.
- Cache repeated genotype patterns within each locus to reduce conversion overhead.
- Prefer optional `python-isal`, then `python-zlib-ng`, for ordinary gzip.
- Use `bgzip -@` for multithreaded BGZF streaming when HTSlib is available.
- Treat `-t` as the total CPU budget and split it between bgzip and parser workers.
- Use TBI/CSI + tabix for true region-parallel BGZF decompression and parsing.
- Preserve indexed contig/window order and filter by VCF POS to avoid overlap duplicates.
- Fall back automatically from a failing indexed mode to BGZF streaming in auto mode.
- Add backend override and troubleshooting options.
- Preserve all original output formats, filtering rules, ordering, and filenames.
