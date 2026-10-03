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
