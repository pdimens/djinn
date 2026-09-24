package fastq

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func writeTempFastq(t *testing.T, name, content string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("writing temp fastq %s: %v", path, err)
	}
	return path
}

func haplotagFastqRecord(id string) string {
	return "@" + id + " VX:i:1\tBX:Z:A01C02B03D04\nACGTACGTAC\n+\nIIIIIIIIII\n"
}

func stlfrFastqRecord(id string) string {
	return "@" + id + "#123_456_789\nACGTACGTAC\n+\nIIIIIIIIII\n"
}

func tellseqFastqRecord(id string) string {
	return "@" + id + ":ATCGATCGATCGATCGAT\nACGTACGTAC\n+\nIIIIIIIIII\n"
}

func TestCheckFastqFormat_Haplotagging(t *testing.T) {
	var sb strings.Builder
	for i := 0; i < 5; i++ {
		sb.WriteString(haplotagFastqRecord("read"))
	}
	path := writeTempFastq(t, "hap.fastq", sb.String())

	fn, err := CheckFastqFormat(path)
	if err != nil {
		t.Fatalf("CheckFastqFormat: %v", err)
	}
	if got, want := string(MissingBarcode), "VX:i:0\tBX:Z:A00C00B00D00"; got != want {
		t.Errorf("MissingBarcode = %q, want %q", got, want)
	}

	reader, err := fastx.NewReader(seq.DNA, path, "")
	if err != nil {
		t.Fatalf("opening reader: %v", err)
	}
	defer reader.Close()
	rec, err := reader.Read()
	if err != nil {
		t.Fatalf("reading record: %v", err)
	}
	bc, valid := fn(rec)
	if bc != "A01C02B03D04" || !valid {
		t.Errorf("fn(rec) = (%q, %v), want (\"A01C02B03D04\", true)", bc, valid)
	}
}

func TestCheckFastqFormat_Stlfr(t *testing.T) {
	var sb strings.Builder
	for i := 0; i < 5; i++ {
		sb.WriteString(stlfrFastqRecord("read"))
	}
	path := writeTempFastq(t, "stlfr.fastq", sb.String())

	fn, err := CheckFastqFormat(path)
	if err != nil {
		t.Fatalf("CheckFastqFormat: %v", err)
	}
	if got, want := string(MissingBarcode), "VX:i:0\tBX:Z:0_0_0"; got != want {
		t.Errorf("MissingBarcode = %q, want %q", got, want)
	}

	reader, err := fastx.NewReader(seq.DNA, path, "")
	if err != nil {
		t.Fatalf("opening reader: %v", err)
	}
	defer reader.Close()
	rec, err := reader.Read()
	if err != nil {
		t.Fatalf("reading record: %v", err)
	}
	bc, valid := fn(rec)
	if bc != "123_456_789" || !valid {
		t.Errorf("fn(rec) = (%q, %v), want (\"123_456_789\", true)", bc, valid)
	}
}

func TestCheckFastqFormat_Tellseq(t *testing.T) {
	var sb strings.Builder
	for i := 0; i < 5; i++ {
		sb.WriteString(tellseqFastqRecord("read"))
	}
	path := writeTempFastq(t, "tellseq.fastq", sb.String())

	fn, err := CheckFastqFormat(path)
	if err != nil {
		t.Fatalf("CheckFastqFormat: %v", err)
	}
	if got, want := string(MissingBarcode), "VX:i:0\tBX:Z:NNNNNNNNNNNNNNNNNN"; got != want {
		t.Errorf("MissingBarcode = %q, want %q", got, want)
	}

	reader, err := fastx.NewReader(seq.DNA, path, "")
	if err != nil {
		t.Fatalf("opening reader: %v", err)
	}
	defer reader.Close()
	rec, err := reader.Read()
	if err != nil {
		t.Fatalf("reading record: %v", err)
	}
	bc, valid := fn(rec)
	if bc != "ATCGATCGATCGATCGAT" || !valid {
		t.Errorf("fn(rec) = (%q, %v), want (\"ATCGATCGATCGATCGAT\", true)", bc, valid)
	}
}

// TestCheckFastqFormat_FewerThan100Records exercises the fix for the bug
// where hitting io.EOF before 100 records were read caused CheckFastqFormat
// to fail outright instead of using the records that were available.
func TestCheckFastqFormat_FewerThan100Records(t *testing.T) {
	var sb strings.Builder
	for i := 0; i < 3; i++ {
		sb.WriteString(stlfrFastqRecord("read"))
	}
	path := writeTempFastq(t, "short.fastq", sb.String())

	fn, err := CheckFastqFormat(path)
	if err != nil {
		t.Fatalf("CheckFastqFormat on a %d-record file returned an error: %v", 3, err)
	}
	if fn == nil {
		t.Fatalf("expected a non-nil detector function")
	}
}

func TestCheckFastqFormat_EmptyFile(t *testing.T) {
	path := writeTempFastq(t, "empty.fastq", "")
	_, err := CheckFastqFormat(path)
	if err == nil {
		t.Fatalf("expected an error for a file with zero records")
	}
}

// TestCheckFastqFormat_MixedFormats exercises the fix for the bug where the
// "more than one linked-read technology format identified" check compared
// (h+s+t) against totalReads, which never fires when different reads (not
// the same read) match different technologies, since the counts still sum
// to totalReads exactly.
func TestCheckFastqFormat_MixedFormats(t *testing.T) {
	var sb strings.Builder
	for i := 0; i < 3; i++ {
		sb.WriteString(haplotagFastqRecord("read"))
	}
	for i := 0; i < 3; i++ {
		sb.WriteString(stlfrFastqRecord("read"))
	}
	path := writeTempFastq(t, "mixed.fastq", sb.String())

	_, err := CheckFastqFormat(path)
	if err == nil {
		t.Fatalf("expected an error when the file mixes haplotagging and stLFR records")
	}
	if !strings.Contains(err.Error(), "more than one linked-read technology") {
		t.Errorf("error = %v, want a message about mixed technology formats", err)
	}
}

func TestCheckFastqFormat_UnrecognizedFormat(t *testing.T) {
	var sb strings.Builder
	for i := 0; i < 5; i++ {
		sb.WriteString("@read\nACGTACGTAC\n+\nIIIIIIIIII\n")
	}
	path := writeTempFastq(t, "unrecognized.fastq", sb.String())

	_, err := CheckFastqFormat(path)
	if err == nil {
		t.Fatalf("expected an error when no linked-read technology can be identified")
	}
}

func TestCheckFastqFormat_MissingFile(t *testing.T) {
	_, err := CheckFastqFormat(filepath.Join(t.TempDir(), "does-not-exist.fastq"))
	if err == nil {
		t.Fatalf("expected an error for a missing file")
	}
}

func TestFqHaplotagBX(t *testing.T) {
	rec := newRecord("read1", "VX:i:1\tBX:Z:A01C02B03D04")
	bc, valid := FqHaplotagBX(rec)
	if bc != "A01C02B03D04" || !valid {
		t.Errorf("FqHaplotagBX = (%q, %v), want (\"A01C02B03D04\", true)", bc, valid)
	}
	// rec must be left untouched (unlike HaplotagBX, this is a read-only probe)
	if got, want := string(rec.Desc), "VX:i:1\tBX:Z:A01C02B03D04"; got != want {
		t.Errorf("rec.Desc was mutated: got %q, want %q", got, want)
	}
}

func TestFqHaplotagBX_NoTag(t *testing.T) {
	rec := newRecord("read1", "no tag here")
	bc, valid := FqHaplotagBX(rec)
	if bc != "" || valid {
		t.Errorf("FqHaplotagBX = (%q, %v), want (\"\", false)", bc, valid)
	}
}

func TestFqStlfr(t *testing.T) {
	rec := newRecord("read1#123_456_789", "")
	bc, valid := FqStlfr(rec)
	if bc != "123_456_789" || !valid {
		t.Errorf("FqStlfr = (%q, %v), want (\"123_456_789\", true)", bc, valid)
	}
	if got, want := string(rec.ID), "read1#123_456_789"; got != want {
		t.Errorf("rec.ID was mutated: got %q, want %q", got, want)
	}
}

func TestFqStlfr_InvalidBarcode(t *testing.T) {
	rec := newRecord("read1#0_456_789", "")
	bc, valid := FqStlfr(rec)
	if bc != "0_456_789" || valid {
		t.Errorf("FqStlfr = (%q, %v), want (\"0_456_789\", false)", bc, valid)
	}
}

func TestFqTellseq(t *testing.T) {
	rec := newRecord("read1:ATCGATCGATCGATCGAT", "")
	bc, valid := FqTellseq(rec)
	if bc != "ATCGATCGATCGATCGAT" || !valid {
		t.Errorf("FqTellseq = (%q, %v), want (\"ATCGATCGATCGATCGAT\", true)", bc, valid)
	}
}
