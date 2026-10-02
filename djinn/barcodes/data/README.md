# 10x_whitelist.txt.zst

The real 10X Genomics barcode whitelist: one 16bp barcode per line,
zstd-compressed, embedded directly into the `djinn` binary via
`//go:embed` (see `../tenxlist.go`). `TenXList`/`NewTenXList()` decodes
and scans it one line at a time -- it never materializes the whole
decompressed list in memory.

The file currently committed here is a tiny synthetic placeholder (a
couple hundred random 16-mers) so the package builds without the real
~4M-line whitelist in this environment. **Replace it with the real file
before shipping** -- same filename, same path, nothing else changes.

## Regenerating it from a real whitelist

`zstdpack.go` in this directory is a standalone helper (excluded from the
normal build via `//go:build ignore`) that recompresses a barcode list
into this format. From the `djinn/djinn` module directory:

```sh
# from a gzip'd source
go run barcodes/data/zstdpack.go -in /path/to/3M-february-2018.txt.gz -out barcodes/data/10x_whitelist.txt.zst

# from a plain text source
go run barcodes/data/zstdpack.go -in /path/to/whitelist.txt -out barcodes/data/10x_whitelist.txt.zst
```

It uses `zstd.WithEncoderLevel(zstd.SpeedBestCompression)` -- slow to
encode (this runs once, offline), fast to decode (this matters, since
`TenXList` scans the whole thing once per `convert` invocation).

## Why zstd over gzip -9

A 16bp barcode over a 4-letter alphabet carries exactly 2 bits/base of
real entropy (32 bits = 4 bytes per barcode), against gzip -9's observed
~16MB for ~4M barcodes (~68MB raw text) -- a ~4.25x ratio. That's already
close to gzip's practical ceiling on this kind of data (DEFLATE's Huffman
coding gets near 2 bits/symbol for a 4-symbol alphabet even though the
input is 8-bit ASCII), so there isn't much more a smarter *general*
compressor can find in the way of long-range redundancy -- whitelists
like this are deliberately chosen to be close to pairwise-random to
maximize Hamming distance between barcodes, which looks like noise to
any compressor.

`zstd` was picked over gzip anyway because:
- It's already a transitive dependency (`klauspost/compress`, pulled in
  via `pgzip`) -- no new dependency added.
- It typically still beats gzip -9 by a meaningful margin on this kind
  of data (bigger window, better entropy coding), even without changing
  the file format.
- Its decoder is significantly faster than gzip's, which matters here
  since decoding happens once per `convert` run over the whole file, not
  once ever.

`xz`/LZMA2 (`ulikunitz/xz`, also already a transitive dependency) can
sometimes squeeze a few more percent out of this exact kind of data than
zstd at max settings, but its decoder is noticeably slower -- not a clearly
better trade here. If shaving every last byte off the shipped binary ever
matters more than decode speed, 2-bit-packing the sequence directly
(A/C/G/T -> 2 bits, no text/newline overhead at all) gets you within a few
percent of the information-theoretic floor (~4 bytes/barcode, ~16MB raw
for 4M barcodes) regardless of compressor -- but that's a custom binary
format needing its own encode/decode code, for a gain that's unlikely to
be more than ~5-10% over zstd on an already near-floor dataset. Worth it
only if binary size is a hard constraint; try `zstd --ultra -22` on the
real file first and see if it's even necessary.
