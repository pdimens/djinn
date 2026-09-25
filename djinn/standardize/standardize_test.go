package standardize

import (
	"io"
	"os"
	"testing"

	"djinn/xam"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func newHeader(t *testing.T) *sam.Header {
	t.Helper()
	hdr, err := sam.NewHeader(nil, nil)
	if err != nil {
		t.Fatalf("sam.NewHeader: %v", err)
	}
	return hdr
}

func newRecord(t *testing.T, name string, bx string, id string) *sam.Record {
	t.Helper()
	rec, err := sam.NewRecord(name, nil, nil, -1, -1, 0, 0, nil, []byte("ACGT"), []byte{40, 40, 40, 40}, nil)
	if err != nil {
		t.Fatalf("sam.NewRecord: %v", err)
	}
	rec.Flags = sam.Unmapped
	if id != "" {
		rec.Name = id
	}
	if bx != "" {
		aux, err := sam.NewAux(sam.Tag{'B', 'X'}, bx)
		if err != nil {
			t.Fatalf("sam.NewAux: %v", err)
		}
		rec.AuxFields = append(rec.AuxFields, aux)
	}
	return rec
}

func writeBamFile(t *testing.T, path string, hdr *sam.Header, recs []*sam.Record) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("os.Create: %v", err)
	}
	defer f.Close()
	bw, err := bam.NewWriterLevel(f, hdr, 0, 1)
	if err != nil {
		t.Fatalf("bam.NewWriterLevel: %v", err)
	}
	for _, r := range recs {
		if err := bw.Write(r); err != nil {
			t.Fatalf("bw.Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
}

func readBamRecords(t *testing.T, path string) []*sam.Record {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("os.Open: %v", err)
	}
	defer f.Close()
	br, err := bam.NewReader(f, 0)
	if err != nil {
		t.Fatalf("bam.NewReader: %v", err)
	}
	var out []*sam.Record
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("br.Read: %v", err)
		}
		out = append(out, rec)
	}
	return out
}

// captureStdout redirects os.Stdout for the duration of fn and returns
// everything written to it.
func captureStdout(t *testing.T, fn func()) []byte {
	t.Helper()
	orig := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("os.Pipe: %v", err)
	}
	os.Stdout = w
	defer func() { os.Stdout = orig }()

	done := make(chan []byte)
	go func() {
		data, _ := io.ReadAll(r)
		done <- data
	}()

	fn()

	w.Close()
	data := <-done
	return data
}

func TestStandardizeAddsBXAndVXFromReadName(t *testing.T) {
	dir := t.TempDir()
	hdr := newHeader(t)

	// stlfr-style barcode embedded in the read name, no existing BX tag.
	recs := []*sam.Record{
		newRecord(t, "read1#1_2_3", "", ""),
	}
	inPath := dir + "/in.bam"
	writeBamFile(t, inPath, hdr, recs)
	outPath := dir + "/out.bam"

	data := captureStdout(t, func() {
		if err := Standardize(inPath, 2, false); err != nil {
			t.Fatalf("Standardize: %v", err)
		}
	})
	if err := os.WriteFile(outPath, data, 0o644); err != nil {
		t.Fatalf("os.WriteFile: %v", err)
	}

	got := readBamRecords(t, outPath)
	if len(got) != 1 {
		t.Fatalf("got %d records, want 1", len(got))
	}
	bx, hasBX := xam.GetStringTag(got[0], "BX")
	if !hasBX {
		t.Fatalf("expected BX tag to be set")
	}
	if bx != "1_2_3" {
		t.Errorf("BX = %q, want %q", bx, "1_2_3")
	}
	vx, hasVX := xam.GetVX(got[0])
	if !hasVX {
		t.Fatalf("expected VX tag to be set")
	}
	if !vx {
		t.Errorf("VX = false, want true for a valid stlfr barcode")
	}
}

func TestStandardizeLeavesRecordsWithoutBarcodeUntouched(t *testing.T) {
	dir := t.TempDir()
	hdr := newHeader(t)

	// No BX tag and a name that doesn't match any recognized barcode pattern.
	recs := []*sam.Record{
		newRecord(t, "plain_read_name", "", ""),
	}
	inPath := dir + "/in.bam"
	writeBamFile(t, inPath, hdr, recs)
	outPath := dir + "/out.bam"

	data := captureStdout(t, func() {
		if err := Standardize(inPath, 2, false); err != nil {
			t.Fatalf("Standardize: %v", err)
		}
	})
	if err := os.WriteFile(outPath, data, 0o644); err != nil {
		t.Fatalf("os.WriteFile: %v", err)
	}

	got := readBamRecords(t, outPath)
	if len(got) != 1 {
		t.Fatalf("got %d records, want 1", len(got))
	}
	if _, hasBX := xam.GetStringTag(got[0], "BX"); hasBX {
		t.Errorf("expected no BX tag to be set when no barcode could be found")
	}
}

func TestStandardizeInvalidBarcodeSetsVXFalse(t *testing.T) {
	dir := t.TempDir()
	hdr := newHeader(t)

	// stlfr barcode containing a 0 segment => invalid per xam.IsValid.
	recs := []*sam.Record{
		newRecord(t, "read1#0_2_3", "", ""),
	}
	inPath := dir + "/in.bam"
	writeBamFile(t, inPath, hdr, recs)
	outPath := dir + "/out.bam"

	data := captureStdout(t, func() {
		if err := Standardize(inPath, 2, false); err != nil {
			t.Fatalf("Standardize: %v", err)
		}
	})
	if err := os.WriteFile(outPath, data, 0o644); err != nil {
		t.Fatalf("os.WriteFile: %v", err)
	}

	got := readBamRecords(t, outPath)
	if len(got) != 1 {
		t.Fatalf("got %d records, want 1", len(got))
	}
	vx, hasVX := xam.GetVX(got[0])
	if !hasVX {
		t.Fatalf("expected VX tag to be set")
	}
	if vx {
		t.Errorf("VX = true, want false for an invalid barcode (leading 0 segment)")
	}
}

func TestStandardizeEmptyInput(t *testing.T) {
	dir := t.TempDir()
	hdr := newHeader(t)

	inPath := dir + "/empty.bam"
	writeBamFile(t, inPath, hdr, nil)
	outPath := dir + "/out.bam"

	data := captureStdout(t, func() {
		if err := Standardize(inPath, 2, false); err != nil {
			t.Fatalf("Standardize on empty input should not error: %v", err)
		}
	})
	if err := os.WriteFile(outPath, data, 0o644); err != nil {
		t.Fatalf("os.WriteFile: %v", err)
	}

	got := readBamRecords(t, outPath)
	if len(got) != 0 {
		t.Errorf("expected 0 records, got %d", len(got))
	}
}

func TestStandardizeLowThreadsClamped(t *testing.T) {
	// threads <= 2 exercises the readThread/writeThread=1/1 branch instead
	// of the >2 split; just verify it still runs correctly end to end.
	dir := t.TempDir()
	hdr := newHeader(t)
	recs := []*sam.Record{newRecord(t, "read1#1_2_3", "", "")}
	inPath := dir + "/in.bam"
	writeBamFile(t, inPath, hdr, recs)
	outPath := dir + "/out.bam"

	data := captureStdout(t, func() {
		if err := Standardize(inPath, 1, false); err != nil {
			t.Fatalf("Standardize: %v", err)
		}
	})
	if err := os.WriteFile(outPath, data, 0o644); err != nil {
		t.Fatalf("os.WriteFile: %v", err)
	}
	got := readBamRecords(t, outPath)
	if len(got) != 1 {
		t.Fatalf("got %d records, want 1", len(got))
	}
}
