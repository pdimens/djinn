package sort

import (
	"container/heap"
	"djinn/xam"
	"errors"
	"io"
	"os"
	"sort"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

const chunkSize = 500_000 // records per temp file, tune to RAM

// SortByBX streams inPath, spills BX:Z-sorted chunks in parallel, k-way merges into outPath.
func SortByBX(infile, outPath, tmpDir string, threads int, asSam bool) error {
	// ── open reader -------------───────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, max(2, threads))

	// ── establish temp dir -------------───────────────────────────────────────────────
	if tmpDir != "" {
		if err := os.MkdirAll(tmpDir, 0o755); err != nil {
			return err
		}
	} else {
		tmpDir = "."
	}

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam sort "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	chunks := make(chan []*sam.Record, threads)
	results := make(chan string, threads)
	errCh := make(chan error, threads)

	var wg sync.WaitGroup
	for range threads {
		wg.Go(func() {
			for c := range chunks {
				path, err := sortAndSpill(c, hdr, tmpDir)
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

	// collector goroutine: drain results as they arrive so workers never block
	var tmpFiles []string
	var collectWG sync.WaitGroup
	collectWG.Go(func() {
		for path := range results {
			tmpFiles = append(tmpFiles, path)
		}
	})

	// producer: read + batch on the calling goroutine
	var chunk []*sam.Record

	for rec := range recChan {
		chunk = append(chunk, rec)
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

	// cleanup temp files on any exit path
	defer func() {
		for _, f := range tmpFiles {
			os.Remove(f)
		}
		os.Remove(tmpDir)
	}()

	select {
	case err := <-errCh:
		return err
	default:
	}

	return mergeChunks(tmpFiles, hdr, outPath, asSam)
}

// bxOf extracts BX:Z; missing tag sorts last (samtools convention: unset tags trail).
func bxOf(r *sam.Record) string {
	if aux, ok := r.Tag([]byte("BX")); ok {
		return aux.Value().(string)
	}
	return "\xff\xff\xff\xff"
}

// sortAndSpill sorts one chunk by BX:Z and writes it to a fresh temp BAM. Each
// call owns its own slice, so no lock is needed across worker goroutines.
func sortAndSpill(chunk []*sam.Record, hdr *sam.Header, tmpDir string) (string, error) {
	sort.Slice(chunk, func(i, j int) bool { return bxOf(chunk[i]) < bxOf(chunk[j]) })

	f, err := os.CreateTemp(tmpDir, ".bxsort-*.bam")
	if err != nil {
		return "", err
	}
	defer f.Close()

	bw, err := bam.NewWriterLevel(f, hdr, 0, 0)
	if err != nil {
		return "", err
	}
	for _, r := range chunk {
		if err := bw.Write(r); err != nil {
			bw.Close()
			return "", err
		}
	}
	if err := bw.Close(); err != nil {
		return "", err
	}
	return f.Name(), nil
}

type mergeItem struct {
	rec *sam.Record
	src int
}
type mergeHeap []mergeItem

func (h mergeHeap) Len() int           { return len(h) }
func (h mergeHeap) Less(i, j int) bool { return bxOf(h[i].rec) < bxOf(h[j].rec) }
func (h mergeHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *mergeHeap) Push(x any)        { *h = append(*h, x.(mergeItem)) }
func (h *mergeHeap) Pop() any {
	old := *h
	n := len(old)
	it := old[n-1]
	*h = old[:n-1]
	return it
}

// mergeChunks does a single-threaded k-way merge of already-sorted temp
// files. This stage is I/O-bound with one writer, so it stays serial;
// parallelizing it would need a second merge pass over partial merges,
// more complexity than the win is worth.
func mergeChunks(tmpFiles []string, hdr *sam.Header, outPath string, asSam bool) (err error) {
	readers := make([]*bam.Reader, len(tmpFiles))
	files := make([]*os.File, len(tmpFiles))
	for i, path := range tmpFiles {
		f, ferr := os.Open(path)
		if ferr != nil {
			return ferr
		}
		files[i] = f
		br, berr := bam.NewReader(f, 0)
		if berr != nil {
			return berr
		}
		readers[i] = br
	}
	defer func() {
		for _, f := range files {
			f.Close()
		}
	}()

	if outPath != "-" {
		out, err := os.Create(outPath)
		if err != nil {
			return err
		}
		defer func() {
			if cerr := out.Close(); err == nil {
				err = cerr
			}
		}()
	}
	// --- Create writer channel --------------------------
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, 2, asSam)

	h := &mergeHeap{}
	for i, br := range readers {
		rec, rerr := br.Read()
		if rerr == nil {
			*h = append(*h, mergeItem{rec, i})
		} else if !errors.Is(rerr, io.EOF) {
			return rerr
		}
	}
	heap.Init(h)

	for h.Len() > 0 {
		top := heap.Pop(h).(mergeItem)
		writeChan <- top.rec
		rec, rerr := readers[top.src].Read()
		if rerr == nil {
			heap.Push(h, mergeItem{rec, top.src})
		} else if !errors.Is(rerr, io.EOF) {
			return rerr
		}
	}
	close(writeChan)
	<-writeDone

	return nil
}
