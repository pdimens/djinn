package convert

import (
	"bufio"
	"bytes"
	"djinn/barcodes"
	"djinn/fastq"
	"fmt"
	"os"
	"strconv"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

// ConvertFqParallel is a parallel-processing variant of ConvertFq. It fans
// each input file's records out to `threads` worker goroutines via
// fastx.Reader.ChunkChan (which Clone()s every record before it crosses
// the channel, so each worker's chunk is independent memory -- see the
// CoreFq doc comment on why that matters), and funnels every worker's
// output through two pieces of explicitly shared, mutex-protected state:
//
//   - seen (which input barcode maps to which already-assigned output
//     barcode) and the barcode Generator itself (which NextInto advances)
//     are inseparable from each other: both are touched only inside
//     assignBarcode, under a single mutex (bcMu), for its whole body.
//     Splitting "check seen" and "call NextInto" across two locks (or two
//     critical sections) would let two goroutines both see a given input
//     barcode as new and hand it two *different* output barcodes --
//     exactly the bug a mutex here exists to prevent.
//   - each output file's *bufio.Writer is shared by every worker
//     processing that file, protected by its own mutex (writeMu). To keep
//     lock contention low, each worker does its CPU-bound work (parsing +
//     barcode assignment + formatting) for an *entire chunk* into a local,
//     unshared buffer first, and only holds writeMu for the single Write
//     of that already-formatted chunk -- one lock/unlock per chunk, not
//     per record.
//
// Caveat -- output ordering: chunks are written in whatever order their
// worker finishes, not necessarily the order they were read in. For a
// single file that's harmless (every record's content is still correct,
// just possibly reordered within the file). For paired R1/R2 conversion,
// though, this file processes R1 and R2 independently and in parallel, so
// a record's position in the output R1 file is not guaranteed to line up
// with its mate's position in the output R2 file the way ConvertFq's
// strictly sequential, one-record-at-a-time loop guarantees. Do not use
// this for a technology/pipeline stage that relies on positional R1/R2
// pairing without adding a chunk-ID-ordered reassembly stage first.
//
// This is an exploratory parallel implementation kept alongside, not in
// place of, ConvertFq.
func ConvertFqParallel(fqs []string, convTo, prefix, bcmap string, threads int) error {
	var err error
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bufio.Writer)

	// --- assess the data type and get the parser -----------------------------
	coreParser, fmtType, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		if bcmap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode file must be provided to identify inline barcodes.", fqs[0])
		}
		return fmt.Errorf("parallel 10X conversion is not implemented in this exploratory version; use ConvertFrom10X")
	}
	if fmtType == convTo {
		return fmt.Errorf("trying to convert to and from identical format (%s)", convTo)
	}

	// --- create the barcode conversion generator -------------------------------------
	switch convTo {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
		converter = (*fastq.CoreFq).ToHaplotagging
	case "stlfr":
		bcs = barcodes.NewStlfr()
		converter = (*fastq.CoreFq).ToStlfr
	case "tellseq":
		bcs = barcodes.NewTellseq()
		converter = (*fastq.CoreFq).ToTellseq
	case "10x":
		bcs, err = barcodes.NewTenXList()
		if err != nil {
			return err
		}
		converter = (*fastq.CoreFq).ToTenX
	default:
		return fmt.Errorf("unknown barcode type %q", convTo)
	}
	defer bcs.Close()

	// ---- shared state, guarded by bcMu -- see the doc comment above ----
	var bcMu sync.Mutex
	seen := make(map[string][]byte, 4_000_000)
	bcBuf := make([]byte, bcs.MaxLen())

	assignBarcode := func(inputBC []byte) ([]byte, error) {
		bcMu.Lock()
		defer bcMu.Unlock()
		if out, ok := seen[string(inputBC)]; ok {
			return out, nil
		}
		n, ok := bcs.NextInto(bcBuf)
		if !ok {
			return nil, fmt.Errorf("too many unique barcodes for the conversion technology requested — unable to generate more barcodes.")
		}
		// independent copy; bcBuf gets reused/overwritten on the next
		// NextInto call, including ones from other goroutines
		out := append([]byte(nil), bcBuf[:n]...)
		seen[string(inputBC)] = out
		return out, nil
	}

	workers := threads
	if workers < 1 {
		workers = 1
	}

	if len(fqs) == 2 {
		// Paired R1/R2: process both files through a chunk-pairing
		// dispatcher so output stays positionally aligned between them --
		// see convertPairedFilesParallel's doc comment.
		r1Out := prefix + ".R1.fq.gz"
		r2Out := prefix + ".R2.fq.gz"
		if err := convertPairedFilesParallel(fqs[0], fqs[1], r1Out, r2Out, workers, coreParser, bcs, converter, assignBarcode); err != nil {
			return err
		}
	} else {
		for idx, fqPath := range fqs {
			outPath := prefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz"
			if err := convertFileParallel(fqPath, outPath, workers, coreParser, bcs, converter, assignBarcode); err != nil {
				return err
			}
		}
	}

	// ------ write barcodes to map ------------------------
	f, err := os.Create(prefix + ".bc.map")
	if err != nil {
		return err
	}
	defer f.Close()

	writer := bufio.NewWriter(f)
	defer writer.Flush()

	for key, val := range seen {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.Write(val)
		writer.WriteByte('\n')
	}

	return nil
}

// chunkSize is the number of records processed as one unit between writer
// lock acquisitions. Larger chunks mean fewer lock/unlock cycles (less
// contention on writeMu) but coarser-grained parallelism and a bigger
// per-chunk memory footprint; this is a reasonable middle ground, not a
// tuned constant.
const convertFqChunkSize = 2000

// convertFileParallel streams fqPath through `workers` goroutines and
// writes the converted output to outPath.
func convertFileParallel(
	fqPath, outPath string,
	workers int,
	coreParser func(*fastx.Record) (fastq.CoreFq, bool),
	bcs barcodes.Generator,
	converter func(*fastq.CoreFq, *bufio.Writer),
	assignBarcode func([]byte) ([]byte, error),
) error {
	fqReader, err := fastx.NewReader(seq.DNA, fqPath, "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqPath, err)
	}
	defer fqReader.Close()

	outfq, err := xopen.Wopen(outPath)
	if err != nil {
		return err
	}
	defer outfq.Close()

	// Buffered enough that the (single) reader goroutine inside ChunkChan
	// can stay ahead of the workers without blocking on every send.
	chunks := fqReader.ChunkChan(workers*2, convertFqChunkSize)

	var writeMu sync.Mutex
	var wg sync.WaitGroup
	errCh := make(chan error, workers)

	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chunk := range chunks {
				if chunk.Err != nil {
					errCh <- fmt.Errorf("reading %s: %w", fqPath, chunk.Err)
					return
				}
				if err := processChunk(chunk.Data, coreParser, bcs, converter, assignBarcode, outfq.Writer, &writeMu); err != nil {
					errCh <- err
					return
				}
			}
		}()
	}

	wg.Wait()
	close(errCh)
	for e := range errCh {
		if e != nil {
			return e
		}
	}
	return nil
}

// formatChunk does the CPU-bound work (barcode parsing/assignment and
// output formatting) for a whole chunk into a local, unshared buffer, and
// returns it unwritten. Touches no shared state except assignBarcode's own
// internal locking -- safe to call from any number of goroutines at once.
func formatChunk(
	records []*fastx.Record,
	coreParser func(*fastx.Record) (fastq.CoreFq, bool),
	bcs barcodes.Generator,
	converter func(*fastq.CoreFq, *bufio.Writer),
	assignBarcode func([]byte) ([]byte, error),
) (*bytes.Buffer, error) {
	var local bytes.Buffer
	bw := bufio.NewWriter(&local)

	for _, rec := range records {
		core, ok := coreParser(rec)
		if !ok {
			// GetInvalid() returns a slice of the generator's precomputed,
			// immutable sentinel value -- safe to call without bcMu,
			// unlike NextInto, which advances mutable generator state.
			core.Barcode = bcs.GetInvalid()
			converter(&core, bw)
			continue
		}
		newBC, err := assignBarcode(core.Barcode)
		if err != nil {
			return nil, err
		}
		core.Barcode = newBC
		converter(&core, bw)
	}
	if err := bw.Flush(); err != nil {
		return nil, err
	}
	return &local, nil
}

// processChunk formats a single-file chunk and writes it to w, taking
// writeMu just long enough for that one Write -- one lock/unlock per
// chunk rather than per record. Used by convertFileParallel, where there
// is only one output file and so no risk of one file's writes outpacing
// another's (see convertPairedFilesParallel for why the paired case needs
// a different locking shape, not this function).
func processChunk(
	records []*fastx.Record,
	coreParser func(*fastx.Record) (fastq.CoreFq, bool),
	bcs barcodes.Generator,
	converter func(*fastq.CoreFq, *bufio.Writer),
	assignBarcode func([]byte) ([]byte, error),
	w *bufio.Writer,
	writeMu *sync.Mutex,
) error {
	local, err := formatChunk(records, coreParser, bcs, converter, assignBarcode)
	if err != nil {
		return err
	}
	writeMu.Lock()
	defer writeMu.Unlock()
	_, err = w.Write(local.Bytes())
	return err
}

// chunkPair is one matched unit of work: chunk k read from R1 together
// with chunk k read from R2. mismatched is set when the dispatcher
// detects the two files don't have the same number of chunks/records --
// in that case Data on whichever side ran out may be nil/short, and the
// pair exists only to carry the error to a worker.
type chunkPair struct {
	r1, r2     fastx.RecordChunk
	mismatched bool
}

// convertPairedFilesParallel processes R1 and R2 together so that output
// stays positionally aligned between the two files, even though absolute
// output order (relative to input order) is not guaranteed -- see
// ConvertFqParallel's doc comment on why that distinction is safe.
//
// A single dispatcher goroutine is the *only* reader of either file's
// ChunkChan. That's the crux of why this works: ChunkChan emits chunk 0,
// 1, 2, ... in strict read order, so chunk k from R1's channel and chunk
// k from R2's channel cover the same record-index range by construction.
// If multiple goroutines raced to pull "the next chunk" from each channel
// independently, worker A could end up with R1's chunk 3 paired against
// whatever R2 chunk happened to be available when it got around to
// receiving -- not necessarily chunk 3. Routing both channels through one
// dispatcher that always does `c1 := <-ch1; c2 := <-ch2` in that fixed
// order removes that race entirely: there is no other reader to compete
// with. The dispatcher then hands each matched pair, as a unit, to a
// worker pool -- so the parallelism is in processing already-correct
// pairs, not in figuring out which chunks match.
//
// Having the same worker *process* both halves of a pair is not, on its
// own, enough to keep R1/R2 output aligned -- an earlier version of this
// function wrote each half under its own independent mutex (writeMu1,
// writeMu2) and that was still broken: nothing stopped another worker's
// pair from writing to w1 (or w2) in the gap between this worker's w1
// write and its w2 write, which lets pair 5's R1 half land before pair
// 3's R2 half even though pair 3's R1 half landed first -- a cross-file
// misalignment with no data race involved, just two independently-locked
// critical sections that were supposed to stay in lockstep and didn't.
// A single mutex (pairWriteMu) held across *both* writes for a pair
// closes that gap: while one worker holds it, no other worker can write
// to either file, so a pair's R1 half and R2 half are always adjacent, in
// the same relative order, in both output files.
func convertPairedFilesParallel(
	r1Path, r2Path, r1Out, r2Out string,
	workers int,
	coreParser func(*fastx.Record) (fastq.CoreFq, bool),
	bcs barcodes.Generator,
	converter func(*fastq.CoreFq, *bufio.Writer),
	assignBarcode func([]byte) ([]byte, error),
) error {
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

	w1, err := xopen.Wopen(r1Out)
	if err != nil {
		return err
	}
	defer w1.Close()

	w2, err := xopen.Wopen(r2Out)
	if err != nil {
		return err
	}
	defer w2.Close()

	ch1 := r1Reader.ChunkChan(workers*2, convertFqChunkSize)
	ch2 := r2Reader.ChunkChan(workers*2, convertFqChunkSize)

	pairs := make(chan chunkPair, workers*2)
	go func() {
		defer close(pairs)
		for {
			c1, ok1 := <-ch1
			c2, ok2 := <-ch2
			switch {
			case !ok1 && !ok2:
				// both files exhausted at the same chunk count -- done
				return
			case ok1 != ok2:
				// one file ran out of chunks before the other: different
				// record counts between R1 and R2.
				pairs <- chunkPair{r1: c1, r2: c2, mismatched: true}
				return
			case len(c1.Data) != len(c2.Data):
				// same chunk count, but the final (partial) chunk sizes
				// differ -- also a record-count mismatch.
				pairs <- chunkPair{r1: c1, r2: c2, mismatched: true}
				return
			default:
				pairs <- chunkPair{r1: c1, r2: c2}
			}
		}
	}()

	// One mutex covering both files' writes for a pair -- see the doc
	// comment above for why two independent per-file mutexes are not
	// sufficient here.
	var pairWriteMu sync.Mutex
	var wg sync.WaitGroup
	errCh := make(chan error, workers)

	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for pair := range pairs {
				if pair.mismatched {
					errCh <- fmt.Errorf("%s and %s do not have the same number of records; paired FASTQ input must be aligned 1:1", r1Path, r2Path)
					return
				}
				if pair.r1.Err != nil {
					errCh <- fmt.Errorf("reading %s: %w", r1Path, pair.r1.Err)
					return
				}
				if pair.r2.Err != nil {
					errCh <- fmt.Errorf("reading %s: %w", r2Path, pair.r2.Err)
					return
				}

				// Format both halves (CPU-bound, no shared state besides
				// assignBarcode's own locking) before taking the write
				// lock, so the lock is held only for the two Writes.
				buf1, err := formatChunk(pair.r1.Data, coreParser, bcs, converter, assignBarcode)
				if err != nil {
					errCh <- err
					return
				}
				buf2, err := formatChunk(pair.r2.Data, coreParser, bcs, converter, assignBarcode)
				if err != nil {
					errCh <- err
					return
				}

				pairWriteMu.Lock()
				_, err1 := w1.Writer.Write(buf1.Bytes())
				_, err2 := w2.Writer.Write(buf2.Bytes())
				pairWriteMu.Unlock()
				if err1 != nil {
					errCh <- err1
					return
				}
				if err2 != nil {
					errCh <- err2
					return
				}
			}
		}()
	}

	wg.Wait()
	close(errCh)
	for e := range errCh {
		if e != nil {
			return e
		}
	}
	return nil
}
