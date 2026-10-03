package convert

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// writeHaplotaggingFixture writes nRecords haplotagging-format FASTQ
// records to path, cycling through nBarcodes distinct BX:Z: barcodes so
// the output exercises both the "already seen" and "assign a new one"
// paths in assignBarcode, across multiple chunks/goroutines.
func writeHaplotaggingFixture(t *testing.T, path string, nRecords, nBarcodes int) []string {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("creating fixture: %v", err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	barcodes := make([]string, nBarcodes)
	for i := range barcodes {
		barcodes[i] = fmt.Sprintf("A%02dC%02dB%02dD%02d", (i%96)+1, ((i*7)%96)+1, ((i*13)%96)+1, ((i*19)%96)+1)
	}

	for i := 0; i < nRecords; i++ {
		bc := barcodes[i%nBarcodes]
		fmt.Fprintf(w, "@read%06d 1:N:0:ATCG\tBX:Z:%s\n", i, bc)
		fmt.Fprintf(w, "ACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")
	}
	return barcodes
}

// writePairedHaplotaggingFixtures writes matching R1/R2 haplotagging FASTQ
// files: record i in each file shares the same base read name (CASAVA
// marker and BX:Z: tag strip off independently per file), so that after
// conversion, a correct implementation must produce record i's R1 and R2
// outputs with the same (stripped) ID at whatever shared position they
// land in, for every i.
func writePairedHaplotaggingFixtures(t *testing.T, r1Path, r2Path string, nRecords, nBarcodes int) {
	t.Helper()
	barcodes := make([]string, nBarcodes)
	for i := range barcodes {
		barcodes[i] = fmt.Sprintf("A%02dC%02dB%02dD%02d", (i%96)+1, ((i*7)%96)+1, ((i*13)%96)+1, ((i*19)%96)+1)
	}

	write := func(path string, casava byte) {
		f, err := os.Create(path)
		if err != nil {
			t.Fatalf("creating fixture: %v", err)
		}
		defer f.Close()
		w := bufio.NewWriter(f)
		defer w.Flush()
		for i := 0; i < nRecords; i++ {
			bc := barcodes[i%nBarcodes]
			fmt.Fprintf(w, "@pair%06d %c:N:0:ATCG\tBX:Z:%s\n", i, casava, bc)
			fmt.Fprintf(w, "ACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")
		}
	}
	write(r1Path, '1')
	write(r2Path, '2')
}

// readFastqIDs parses a converted tellseq-format FASTQ ("@<id>:<barcode>")
// and returns the base read ID (barcode stripped) of every record, in
// file order.
func readFastqIDs(t *testing.T, path string) []string {
	t.Helper()
	r, err := fastx.NewReader(seq.DNA, path, "")
	if err != nil {
		t.Fatalf("opening %s: %v", path, err)
	}
	defer r.Close()

	var out []string
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		id := string(rec.ID)
		for i := len(id) - 1; i >= 0; i-- {
			if id[i] == ':' {
				out = append(out, id[:i])
				break
			}
		}
	}
	return out
}

// readFastqBarcodes parses a converted tellseq-format FASTQ and returns
// the inline barcode of every record, in file order.
func readFastqBarcodes(t *testing.T, path string) []string {
	t.Helper()
	r, err := fastx.NewReader(seq.DNA, path, "")
	if err != nil {
		t.Fatalf("opening %s: %v", path, err)
	}
	defer r.Close()

	var out []string
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		// tellseq format: "@<id>:<barcode>\t<casava>"
		id := string(rec.ID)
		for i := len(id) - 1; i >= 0; i-- {
			if id[i] == ':' {
				out = append(out, id[i+1:])
				break
			}
		}
	}
	return out
}

func readBcMap(t *testing.T, path string) map[string]string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("reading %s: %v", path, err)
	}
	m := make(map[string]string)
	for _, line := range splitLines(string(data)) {
		if line == "" {
			continue
		}
		var k, v string
		for i := 0; i < len(line); i++ {
			if line[i] == '\t' {
				k, v = line[:i], line[i+1:]
				break
			}
		}
		m[k] = v
	}
	return m
}

func splitLines(s string) []string {
	var out []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == '\n' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	if start < len(s) {
		out = append(out, s[start:])
	}
	return out
}

func TestConvertFqParallel_CorrectnessUnderConcurrency(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fq")
	const nRecords = 6000 // spans 3 chunks at convertFqChunkSize=2000
	const nBarcodes = 50
	writeHaplotaggingFixture(t, in, nRecords, nBarcodes)

	prefix := filepath.Join(dir, "out")
	if err := ConvertFqParallel([]string{in}, "tellseq", prefix, "", 4); err != nil {
		t.Fatalf("ConvertFqParallel: %v", err)
	}

	gotBarcodes := readFastqBarcodes(t, prefix+".R1.fq.gz")
	if len(gotBarcodes) != nRecords {
		t.Fatalf("got %d output records, want %d", len(gotBarcodes), nRecords)
	}

	bcMap := readBcMap(t, prefix+".bc.map")
	if len(bcMap) != nBarcodes {
		t.Fatalf("bc.map has %d entries, want %d (one per distinct input barcode)", len(bcMap), nBarcodes)
	}

	// Every input barcode must have mapped to exactly one output barcode
	// throughout the run, and distinct input barcodes must never collide
	// on the same output barcode.
	seenOutputs := make(map[string]string) // output -> input, to catch collisions
	for input, output := range bcMap {
		if other, collided := seenOutputs[output]; collided && other != input {
			t.Fatalf("output barcode %q assigned to two different input barcodes: %q and %q", output, other, input)
		}
		seenOutputs[output] = input
	}

	// Cross-check: every barcode actually written into the output file
	// must be one of the values in bc.map (i.e. the writer used the same
	// assignment assignBarcode produced, not something stale/racy).
	validOutputs := make(map[string]bool, len(bcMap))
	for _, output := range bcMap {
		validOutputs[output] = true
	}
	for i, bc := range gotBarcodes {
		if !validOutputs[bc] {
			t.Fatalf("record %d has output barcode %q, which is not in bc.map", i, bc)
		}
	}
}

func TestConvertFqParallel_SingleWorkerMatchesSequentialAssignmentCount(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fq")
	const nRecords = 500
	const nBarcodes = 10
	writeHaplotaggingFixture(t, in, nRecords, nBarcodes)

	seqPrefix := filepath.Join(dir, "seq")
	if err := ConvertFq([]string{in}, "tellseq", seqPrefix, "", 1); err != nil {
		t.Fatalf("ConvertFq: %v", err)
	}
	parPrefix := filepath.Join(dir, "par")
	if err := ConvertFqParallel([]string{in}, "tellseq", parPrefix, "", 4); err != nil {
		t.Fatalf("ConvertFqParallel: %v", err)
	}

	seqMap := readBcMap(t, seqPrefix+".bc.map")
	parMap := readBcMap(t, parPrefix+".bc.map")
	if len(seqMap) != nBarcodes || len(parMap) != nBarcodes {
		t.Fatalf("expected %d distinct barcodes in both maps, got sequential=%d parallel=%d", nBarcodes, len(seqMap), len(parMap))
	}
	// Note: *which* output barcode a given input barcode receives can
	// legitimately differ between the two (concurrent "first encounter"
	// order isn't the same as sequential read order) -- only the shape
	// (same count, no collisions, same record count) is guaranteed equal.
	seqBarcodes := readFastqBarcodes(t, seqPrefix+".R1.fq.gz")
	parBarcodes := readFastqBarcodes(t, parPrefix+".R1.fq.gz")
	if len(seqBarcodes) != len(parBarcodes) {
		t.Fatalf("record count differs: sequential=%d parallel=%d", len(seqBarcodes), len(parBarcodes))
	}
}

func TestConvertFqParallel_UnknownFormat(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fq")
	writeHaplotaggingFixture(t, in, 200, 5)

	err := ConvertFqParallel([]string{in}, "not-a-real-format", filepath.Join(dir, "out"), "", 4)
	if err == nil {
		t.Fatal("expected an error for an unknown target format")
	}
}

// This is the test that actually exercises the chunk-pairing dispatcher:
// with enough records to span many chunks and enough workers that chunk
// pairs can plausibly complete out of their original order, R1 and R2
// output must still agree record-for-record at every shared position.
// Without the single-dispatcher design (e.g. if two independent worker
// pools pulled from R1's and R2's ChunkChan separately), this is exactly
// the test that would catch the resulting misalignment.
func TestConvertFqParallel_PairedOutputStaysAligned(t *testing.T) {
	dir := t.TempDir()
	r1In := filepath.Join(dir, "r1.fq")
	r2In := filepath.Join(dir, "r2.fq")
	const nRecords = 12000 // 6 chunks at convertFqChunkSize=2000
	const nBarcodes = 37
	writePairedHaplotaggingFixtures(t, r1In, r2In, nRecords, nBarcodes)

	prefix := filepath.Join(dir, "out")
	if err := ConvertFqParallel([]string{r1In, r2In}, "tellseq", prefix, "", 8); err != nil {
		t.Fatalf("ConvertFqParallel: %v", err)
	}

	r1IDs := readFastqIDs(t, prefix+".R1.fq.gz")
	r2IDs := readFastqIDs(t, prefix+".R2.fq.gz")

	if len(r1IDs) != nRecords || len(r2IDs) != nRecords {
		t.Fatalf("got %d R1 records and %d R2 records, want %d each", len(r1IDs), len(r2IDs), nRecords)
	}
	for i := range r1IDs {
		if r1IDs[i] != r2IDs[i] {
			t.Fatalf("position %d misaligned: R1 ID %q != R2 ID %q", i, r1IDs[i], r2IDs[i])
		}
	}

	// Confirm the fixture actually got reordered relative to input at
	// least somewhere -- otherwise this test wouldn't be exercising
	// anything the single-file case doesn't already cover. (Not a
	// correctness requirement, just a sanity check that this test has
	// teeth: enough chunks/workers that completion order isn't trivially
	// already-sequential.)
	reordered := false
	for i, id := range r1IDs {
		if id != fmt.Sprintf("pair%06d", i) {
			reordered = true
			break
		}
	}
	t.Logf("output reordered relative to input: %v", reordered)
}
