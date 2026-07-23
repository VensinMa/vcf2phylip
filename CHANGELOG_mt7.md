# 2.9-mt7

- Added direct byte-range parallel reading for uncompressed VCF files.
- Each worker independently aligns nominal offsets to complete VCF records.
- Uses `os.pread()` where available and falls back to `seek()`/`read()`.
- Removes parent-process VCF chunk serialization for plain multi-worker input.
- Keeps results ordered by byte range and preserves all existing matrix formats.
- Added `--input-backend plain-stream` to force the mt6-style input path.
- Added fallback from direct range mode to sequential plain streaming.
- Added tests for range boundaries inside a multi-megabyte VCF record.
