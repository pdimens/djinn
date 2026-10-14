// Package fqsort sorts paired-end FASTQ reads by linked-read barcode,
// mirroring the sort package's BAM/SAM k-way merge sort: read + pair,
// sort-and-spill chunks in parallel, then k-way merge the spilled chunks
// back into barcode order.
package fqsort

import (
	"bufio"
	"bytes"
	"djinn/barcodes"
	"djinn/fastq"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"slices"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

const chunkSize = 250_000 // read pairs per temp file, tune to RAM

// missingBXKey sorts after every real barcode (same samtools-derived
// convention as sort.missingBXKey). It's only ever used as a sort key,
// never written out: a pair with no resolvable barcode keeps a nil
// fqPair.barcode, so ToStandard emits the same empty "BX:Z:" it would for
// any other caller that never found one (fastq.Haplotag2Corefq, etc.).
var missingBXKey = []byte{0xff, 0xff, 0xff, 0xff}

// fqPair is a self-contained (deep-copied) read pair: both mates' FASTQ
// fields plus the barcode/validity resolved from whichever mate carried
// it. Unlike the xam sort's bxKeyOf, which can alias the BAM reader's own
// buffers because a record is re-encoded to its temp file immediately,
// fqPair must own its memory: a chunk of these sits in RAM across a sort
// before being spilled, long after fastx.Reader has reused its internal
// buffers for later records.
//
// id and comments are stored once, not once per mate: every format
// CoreFqParser detects carries its barcode/CASAVA marker as the only
// thing that structurally differs between mates, so once those are
// stripped out (which CoreFq already does), the remaining description
// text is identical on both -- see resolveComments. CASAVA itself isn't
// stored at all: ToStandard only ever emits its first byte (the '1'/'2'
// mate digit, see casavaOne/casavaTwo below), so all a stored pair needs
// is whether each mate had a marker at all.
type fqPair struct {
	barcode                []byte // nil if neither mate carried one
	id                     []byte // shared base read ID (CASAVA suffix stripped)
	comments               []byte // shared description text once CASAVA/barcode are stripped
	valid                  bool
	hasCasava1, hasCasava2 bool
	seq1, qual1            []byte
	seq2, qual2            []byte
}

// pairKey returns p's sort key: its barcode, or missingBXKey when p has
// none (unset barcodes sort last, matching the xam sort's convention).
func pairKey(p *fqPair) []byte {
	if len(p.barcode) == 0 {
		return missingBXKey
	}
	return p.barcode
}

// resolveBarcode prefers mate 1's barcode, falling back to mate 2's; nil
// if neither mate's core parsed to one.
func resolveBarcode(c1 fastq.CoreFq, ok1 bool, c2 fastq.CoreFq, ok2 bool) []byte {
	if ok1 && len(c1.Barcode) > 0 {
		return c1.Barcode
	}
	if ok2 && len(c2.Barcode) > 0 {
		return c2.Barcode
	}
	return nil
}

// resolveValid mirrors resolveBarcode's preference, since a record with
// no parsed core has a meaningless Valid value.
func resolveValid(c1 fastq.CoreFq, ok1 bool, c2 fastq.CoreFq, ok2 bool) bool {
	if ok1 {
		return c1.Valid
	}
	if ok2 {
		return c2.Valid
	}
	return false
}

// resolveComments assumes both mates carry identical description text
// once their barcode/CASAVA markers are stripped -- true for every format
// CoreFqParser recognizes, since the barcode tag/marker is the only thing
// that differs structurally between mates. Prefers whichever side is
// non-empty, defaulting to mate 1, so that assumption's one failure mode
// (the two genuinely differ) still keeps *something* rather than
// silently losing both.
func resolveComments(c1, c2 fastq.CoreFq) []byte {
	if len(c1.Comments) > 0 {
		return c1.Comments
	}
	return c2.Comments
}

// casavaOne and casavaTwo are the literal CASAVA bytes ToStandard ever
// writes -- it only ever reads rec.CASAVA[0], the mate digit -- so a
// stored pair only needs to know whether each mate had a marker at all,
// not the marker's own text.
var casavaOne = []byte{'1'}
var casavaTwo = []byte{'2'}

// SortFqByBX reads r1Path/r2Path as a matched read-pair stream, sorts the
// pairs by linked-read barcode, and writes outPrefix+".R1.fq.gz" /
// outPrefix+".R2.fq.gz". The barcode format (haplotagging/stLFR/TELL-seq)
// is auto-detected from r1Path, same as fastq.CoreFqParser's other
// callers. When none of those is detected and bcMap names a 10X barcode
// whitelist file, falls back to the 10X path instead -- same dispatch
// convert.ConvertFq uses: format detection first, an explicit whitelist
// as the fallback, since 10X's barcode lives inline in R1's sequence and
// has no marker a sampling pass could recognize on its own.
//
// Output is written in the same format it was read in (haplotagging stays
// haplotagging, etc.) unless standardize is true, in which case output
// uses fastq.CoreFq.ToStandard's generic BX:Z:/VX:i: form instead -- the
// same opt-in every other djinn fastq tool (count, extract, invalid,
// sample, singletons) already leaves to the caller.
func SortFqByBX(r1Path, r2Path, outPrefix, tmpDir string, threads int, standardize bool, bcMap string) error {
	seq.ValidateSeq = false

	coreParser, fmtType, err := fastq.CoreFqParser(r1Path)
	if err != nil {
		return err
	}

	// parseMate1/parseMate2 are asymmetric only for 10X, where the
	// barcode lives solely in R1's sequence prefix: TenXR22Corefq never
	// finds one on its own, so its ok is always false, which is exactly
	// what resolveBarcode/resolveValid need to fall through to mate 1's
	// result below -- the same preference they already apply for every
	// other format's (rarer, but possible) single-mate-only barcode.
	parseMate1, parseMate2 := coreParser, coreParser
	if coreParser == nil {
		if bcMap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode whitelist file must be provided to identify inline barcodes.", r1Path)
		}
		bclist, err := barcodes.LoadTenXFromFile(bcMap)
		if err != nil {
			return err
		}
		parseMate1 = func(rec *fastx.Record) (fastq.CoreFq, bool) { return fastq.TenXR12Corefq(rec, bclist) }
		parseMate2 = func(rec *fastx.Record) (fastq.CoreFq, bool) { return fastq.TenXR22Corefq(rec), false }
		fmtType = "10x"
	}

	converter, err := formatConverter(fmtType, standardize)
	if err != nil {
		return err
	}

	r1Reader, err := fastx.NewReader(seq.DNA, r1Path, "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", r1Path, err)
	}
	defer r1Reader.Close()

	r2Reader, err := fastx.NewReader(seq.DNA, r2Path, "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", r2Path, err)
	}
	defer r2Reader.Close()

	if tmpDir != "" {
		if err := os.MkdirAll(tmpDir, 0o755); err != nil {
			return err
		}
	} else {
		tmpDir = "."
	}

	chunks := make(chan []fqPair, threads)
	results := make(chan string, threads)
	errCh := make(chan error, threads)

	var wg sync.WaitGroup
	for range threads {
		wg.Go(func() {
			for c := range chunks {
				path, err := sortAndSpill(c, tmpDir)
				if err != nil {
					select {
					case errCh <- err:
					default:
					}
					continue
				}
				results <- path
			}
		})
	}

	var tmpFiles []string
	var collectWG sync.WaitGroup
	collectWG.Go(func() {
		for path := range results {
			tmpFiles = append(tmpFiles, path)
		}
	})

	// producer: read matched pairs on the calling goroutine (reading is
	// inherently serial here -- R1 and R2 must be advanced in lockstep --
	// so there's no analogue to the xam sort's separate reader channel).
	var readErr error
	var chunk []fqPair
readLoop:
	for {
		rec1, err1 := r1Reader.Read()
		rec2, err2 := r2Reader.Read()
		eof1, eof2 := err1 == io.EOF, err2 == io.EOF
		if eof1 && eof2 {
			break readLoop
		}
		if err1 != nil && !eof1 {
			readErr = err1
			break readLoop
		}
		if err2 != nil && !eof2 {
			readErr = err2
			break readLoop
		}
		if eof1 != eof2 {
			readErr = fmt.Errorf("%s and %s have different numbers of records (mate files out of sync)", r1Path, r2Path)
			break readLoop
		}

		core1, ok1 := parseMate1(rec1)
		core2, ok2 := parseMate2(rec2)
		if !bytes.Equal(core1.ID, core2.ID) {
			readErr = fmt.Errorf("mismatched read IDs between mates: %q (R1) vs %q (R2)", core1.ID, core2.ID)
			break readLoop
		}

		chunk = append(chunk, fqPair{
			barcode:    append([]byte(nil), resolveBarcode(core1, ok1, core2, ok2)...),
			id:         append([]byte(nil), core1.ID...),
			comments:   append([]byte(nil), resolveComments(core1, core2)...),
			valid:      resolveValid(core1, ok1, core2, ok2),
			hasCasava1: len(core1.CASAVA) > 0,
			hasCasava2: len(core2.CASAVA) > 0,
			seq1:       append([]byte(nil), core1.Seq...),
			qual1:      append([]byte(nil), core1.Qual...),
			seq2:       append([]byte(nil), core2.Seq...),
			qual2:      append([]byte(nil), core2.Qual...),
		})
		if len(chunk) >= chunkSize {
			chunks <- chunk
			chunk = nil
		}
	}

	if len(chunk) > 0 {
		chunks <- chunk
	}
	close(chunks)

	wg.Wait()
	close(results)
	collectWG.Wait()

	defer func() {
		for _, f := range tmpFiles {
			os.Remove(f)
		}
		os.Remove(tmpDir)
	}()

	if readErr != nil {
		return readErr
	}

	select {
	case err := <-errCh:
		return err
	default:
	}

	return mergeFqChunks(tmpFiles, outPrefix, converter, standardize)
}

// formatConverter picks which fastq.CoreFq writer method reconstructs
// output records: the format-specific one matching fmtType (so output
// round-trips the input format), or ToStandard when standardize is set.
func formatConverter(fmtType string, standardize bool) (func(*fastq.CoreFq, *bufio.Writer), error) {
	if standardize {
		return (*fastq.CoreFq).ToStandard, nil
	}
	switch fmtType {
	case "haplotagging":
		return (*fastq.CoreFq).ToHaplotagging, nil
	case "stlfr":
		return (*fastq.CoreFq).ToStlfr, nil
	case "tellseq":
		return (*fastq.CoreFq).ToTellseq, nil
	case "10x":
		return (*fastq.CoreFq).ToTenX, nil
	default:
		return nil, fmt.Errorf("unsupported linked-read technology %q for format-preserving output; use --standardize", fmtType)
	}
}

// sortAndSpill sorts one chunk by barcode and writes it to a fresh,
// uncompressed temp file in djinn's own binary pair format (see
// writePair/readPair) -- not FASTQ text, since this file is only ever
// read back by mergeFqChunks and re-parsing FASTQ text would mean paying
// the barcode/CASAVA extraction cost a second time for no benefit.
func sortAndSpill(chunk []fqPair, tmpDir string) (string, error) {
	slices.SortFunc(chunk, func(a, b fqPair) int { return bytes.Compare(pairKey(&a), pairKey(&b)) })

	f, err := os.CreateTemp(tmpDir, ".fqsort-*.bin")
	if err != nil {
		return "", err
	}
	defer f.Close()

	w := bufio.NewWriterSize(f, 1<<20)
	for i := range chunk {
		if err := writePair(w, &chunk[i]); err != nil {
			return "", err
		}
	}
	if err := w.Flush(); err != nil {
		return "", err
	}
	return f.Name(), nil
}

// ── binary spill format ─────────────────────────────────────────────────────

func writeBytes(w *bufio.Writer, b []byte) error {
	var lenBuf [4]byte
	binary.LittleEndian.PutUint32(lenBuf[:], uint32(len(b)))
	if _, err := w.Write(lenBuf[:]); err != nil {
		return err
	}
	if len(b) > 0 {
		if _, err := w.Write(b); err != nil {
			return err
		}
	}
	return nil
}

func readBytes(r *bufio.Reader) ([]byte, error) {
	var lenBuf [4]byte
	if _, err := io.ReadFull(r, lenBuf[:]); err != nil {
		return nil, err
	}
	n := binary.LittleEndian.Uint32(lenBuf[:])
	if n == 0 {
		return nil, nil
	}
	b := make([]byte, n)
	if _, err := io.ReadFull(r, b); err != nil {
		return nil, err
	}
	return b, nil
}

// flag bits packed into writePair/readPair's single trailing byte.
const (
	flagValid      = 1 << 0
	flagHasCasava1 = 1 << 1
	flagHasCasava2 = 1 << 2
)

// writePair/readPair's field order must match exactly. Barcode leads
// since it's the sort key -- were this format ever inspected line-by-line
// (e.g. hexdump during debugging) the key would be the first thing in
// each record, same spirit as putting it first in a hypothetical TSV
// encoding.
func writePair(w *bufio.Writer, p *fqPair) error {
	for _, b := range [...][]byte{
		p.barcode, p.id, p.comments,
		p.seq1, p.qual1,
		p.seq2, p.qual2,
	} {
		if err := writeBytes(w, b); err != nil {
			return err
		}
	}
	var flags byte
	if p.valid {
		flags |= flagValid
	}
	if p.hasCasava1 {
		flags |= flagHasCasava1
	}
	if p.hasCasava2 {
		flags |= flagHasCasava2
	}
	return w.WriteByte(flags)
}

// readPair reads one record written by writePair. err is io.EOF (and p is
// nil) only when the stream is exhausted exactly at a record boundary, same
// contract sort.mergeChunks relies on for its BAM readers.
func readPair(r *bufio.Reader) (*fqPair, error) {
	var fields [7][]byte
	for i := range fields {
		b, err := readBytes(r)
		if err != nil {
			return nil, err
		}
		fields[i] = b
	}
	flags, err := r.ReadByte()
	if err != nil {
		return nil, err
	}
	return &fqPair{
		barcode:    fields[0],
		id:         fields[1],
		comments:   fields[2],
		seq1:       fields[3],
		qual1:      fields[4],
		seq2:       fields[5],
		qual2:      fields[6],
		valid:      flags&flagValid != 0,
		hasCasava1: flags&flagHasCasava1 != 0,
		hasCasava2: flags&flagHasCasava2 != 0,
	}, nil
}

// ── k-way merge ──────────────────────────────────────────────────────────────

// fqMergeItem is one entry in the merge heap: a pair and which temp file
// it came from. Hand-rolled (not container/heap) for the same reason as
// sort.mergeHeap: avoid boxing every pair into an interface{} on every
// heap operation.
type fqMergeItem struct {
	pair *fqPair
	src  int
}

type fqMergeHeap []fqMergeItem

func (h *fqMergeHeap) push(it fqMergeItem) {
	*h = append(*h, it)
	i := len(*h) - 1
	for i > 0 {
		parent := (i - 1) / 2
		if bytes.Compare(pairKey((*h)[i].pair), pairKey((*h)[parent].pair)) >= 0 {
			break
		}
		(*h)[i], (*h)[parent] = (*h)[parent], (*h)[i]
		i = parent
	}
}

func (h *fqMergeHeap) pop() fqMergeItem {
	old := *h
	n := len(old) - 1
	top := old[0]
	old[0] = old[n]
	*h = old[:n]
	h.siftDown(0)
	return top
}

func (h *fqMergeHeap) siftDown(i int) {
	n := len(*h)
	for {
		left := 2*i + 1
		if left >= n {
			return
		}
		smallest := left
		if right := left + 1; right < n && bytes.Compare(pairKey((*h)[right].pair), pairKey((*h)[left].pair)) < 0 {
			smallest = right
		}
		if bytes.Compare(pairKey((*h)[smallest].pair), pairKey((*h)[i].pair)) >= 0 {
			return
		}
		(*h)[i], (*h)[smallest] = (*h)[smallest], (*h)[i]
		i = smallest
	}
}

// mergeFqChunks does a single-threaded k-way merge of the already
// barcode-sorted temp files, writing the result as gzipped paired FASTQ
// via converter. standardize controls how an absent CASAVA marker is
// handled -- see writePairAsFastq.
func mergeFqChunks(tmpFiles []string, outPrefix string, converter func(*fastq.CoreFq, *bufio.Writer), standardize bool) (err error) {
	readers := make([]*bufio.Reader, len(tmpFiles))
	files := make([]*os.File, len(tmpFiles))
	for i, path := range tmpFiles {
		f, ferr := os.Open(path)
		if ferr != nil {
			return ferr
		}
		files[i] = f
		readers[i] = bufio.NewReaderSize(f, 1<<20)
	}
	defer func() {
		for _, f := range files {
			f.Close()
		}
	}()

	r1Out, err := xopen.Wopen(outPrefix + ".R1.fq.gz")
	if err != nil {
		return err
	}
	defer r1Out.Close()

	r2Out, err := xopen.Wopen(outPrefix + ".R2.fq.gz")
	if err != nil {
		return err
	}
	defer r2Out.Close()

	h := &fqMergeHeap{}
	for i, r := range readers {
		p, perr := readPair(r)
		if perr == nil {
			h.push(fqMergeItem{pair: p, src: i})
		} else if perr != io.EOF {
			return perr
		}
	}

	for len(*h) > 0 {
		top := h.pop()
		writePairAsFastq(top.pair, r1Out.Writer, r2Out.Writer, converter, standardize)
		p, perr := readPair(readers[top.src])
		if perr == nil {
			h.push(fqMergeItem{pair: p, src: top.src})
		} else if perr != io.EOF {
			return perr
		}
	}

	return nil
}

// writePairAsFastq re-assembles both mates as fastq.CoreFq and writes them
// with converter, reproducing either the original format (converter is
// one of ToHaplotagging/ToStlfr/ToTellseq) or the generic BX:Z:/VX:i: form
// every other djinn fastq tool already produces (converter is ToStandard,
// when standardize is set).
//
// Only ToStandard treats an absent CASAVA marker as optional (it omits
// the "/1"-style suffix entirely, matching every other CoreFq-producing
// caller that never found one); the format-specific writers always
// expect *a* mate digit to embed -- ToHaplotagging in particular indexes
// CASAVA[0] unconditionally and would panic on a nil/empty CASAVA -- so
// when standardize is false the mate digit is always synthesized from
// which side of the pair this is, regardless of whether the source
// record actually carried a marker.
func writePairAsFastq(p *fqPair, w1, w2 *bufio.Writer, converter func(*fastq.CoreFq, *bufio.Writer), standardize bool) {
	casava1, casava2 := casavaOne, casavaTwo
	if standardize {
		casava1, casava2 = nil, nil
		if p.hasCasava1 {
			casava1 = casavaOne
		}
		if p.hasCasava2 {
			casava2 = casavaTwo
		}
	}
	c1 := fastq.CoreFq{ID: p.id, Seq: p.seq1, Qual: p.qual1, CASAVA: casava1, Barcode: p.barcode, Valid: p.valid, Comments: p.comments}
	converter(&c1, w1)
	c2 := fastq.CoreFq{ID: p.id, Seq: p.seq2, Qual: p.qual2, CASAVA: casava2, Barcode: p.barcode, Valid: p.valid, Comments: p.comments}
	converter(&c2, w2)
}
