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

	for idx, fqPath := range fqs {
		outPath := prefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz"
		if err := convertFileParallel(fqPath, outPath, workers, coreParser, bcs, converter, assignBarcode); err != nil {
			return err
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

// processChunk does the CPU-bound work (barcode parsing/assignment and
// output formatting) for a whole chunk into a local, unshared buffer, then
// takes writeMu just long enough to flush that buffer to the shared
// output writer -- one lock/unlock per chunk rather than per record.
func processChunk(
	records []*fastx.Record,
	coreParser func(*fastx.Record) (fastq.CoreFq, bool),
	bcs barcodes.Generator,
	converter func(*fastq.CoreFq, *bufio.Writer),
	assignBarcode func([]byte) ([]byte, error),
	w *bufio.Writer,
	writeMu *sync.Mutex,
) error {
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
			return err
		}
		core.Barcode = newBC
		converter(&core, bw)
	}
	if err := bw.Flush(); err != nil {
		return err
	}

	writeMu.Lock()
	defer writeMu.Unlock()
	_, err := w.Write(local.Bytes())
	return err
}
