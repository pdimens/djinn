package convert

import (
	"bufio"
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// captureStdout redirects os.Stdout for the duration of fn and returns
// everything written to it. ConvertXam always writes its BAM/SAM output to
// "-" (stdout), so this is required to exercise it without polluting the
// test process's real stdout.
func captureStdout(t *testing.T, fn func()) []byte {
	t.Helper()
	orig := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("os.Pipe: %v", err)
	}
	os.Stdout = w

	done := make(chan []byte)
	go func() {
		var buf bytes.Buffer
		io.Copy(&buf, r)
		done <- buf.Bytes()
	}()

	fn()

	w.Close()
	os.Stdout = orig
	return <-done
}

func writeSam(t *testing.T, path, content string) {
	t.Helper()
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

const testSamHeader = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"

func TestConvertXamUnknownType(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	writeSam(t, in, testSamHeader+"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n")

	err := ConvertXam(in, "not-a-real-format", filepath.Join(dir, "map.tsv"), 1, true)
	if err == nil {
		t.Fatal("expected error for unknown barcode type, got nil")
	}
}

func TestConvertXamHaplotagging(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	bcMap := filepath.Join(dir, "map.tsv")

	samText := testSamHeader +
		"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A01C01B01D01\n" +
		"read2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A01C01B01D01\n" + // same barcode as read1: should reuse
		"read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n" + // no barcode at all
		"read4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A00C01B01D01\n" // invalid barcode (A00 segment)
	writeSam(t, in, samText)

	var convErr error
	out := captureStdout(t, func() {
		convErr = ConvertXam(in, "haplotagging", bcMap, 1, true)
	})
	if convErr != nil {
		t.Fatalf("ConvertXam: %v", convErr)
	}

	lines := strings.Split(strings.TrimRight(string(out), "\n"), "\n")
	var recLines []string
	for _, l := range lines {
		if strings.HasPrefix(l, "@") {
			continue
		}
		if l == "" {
			continue
		}
		recLines = append(recLines, l)
	}
	if len(recLines) != 4 {
		t.Fatalf("got %d record lines, want 4:\n%s", len(recLines), out)
	}

	bxOf := func(line string) string {
		for _, field := range strings.Split(line, "\t") {
			if strings.HasPrefix(field, "BX:Z:") {
				return strings.TrimPrefix(field, "BX:Z:")
			}
		}
		return ""
	}

	bx1 := bxOf(recLines[0])
	bx2 := bxOf(recLines[1])
	bx3 := bxOf(recLines[2])
	bx4 := bxOf(recLines[3])

	if bx1 == "" {
		t.Fatal("read1 missing converted BX tag")
	}
	if bx1 != bx2 {
		t.Errorf("read1 and read2 shared the same source barcode but got different converted barcodes: %q vs %q", bx1, bx2)
	}
	if bx3 == "" || bx4 == "" {
		t.Fatal("read3/read4 should have been assigned the sentinel invalid barcode")
	}
	if bx3 != bx4 {
		t.Errorf("read3 and read4 (both invalid) should share the sentinel invalid barcode, got %q vs %q", bx3, bx4)
	}
	if bx1 == bx3 {
		t.Errorf("valid-barcode conversion collided with the invalid sentinel barcode: %q", bx1)
	}

	// the barcode map should contain exactly one entry: the one unique,
	// valid source barcode that was actually converted.
	mapContent, err := os.ReadFile(bcMap)
	if err != nil {
		t.Fatalf("ReadFile bcMap: %v", err)
	}
	scanner := bufio.NewScanner(bytes.NewReader(mapContent))
	var mapLines []string
	for scanner.Scan() {
		mapLines = append(mapLines, scanner.Text())
	}
	if len(mapLines) != 1 {
		t.Fatalf("bcMap has %d lines, want 1:\n%s", len(mapLines), mapContent)
	}
	if !strings.HasPrefix(mapLines[0], "A01C01B01D01\t") {
		t.Errorf("bcMap line = %q, want prefix %q", mapLines[0], "A01C01B01D01\t")
	}
}

func TestConvertXamEmptyInput(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	writeSam(t, in, testSamHeader)
	bcMap := filepath.Join(dir, "map.tsv")

	var convErr error
	out := captureStdout(t, func() {
		convErr = ConvertXam(in, "haplotagging", bcMap, 1, true)
	})
	if convErr != nil {
		t.Fatalf("ConvertXam: %v", convErr)
	}
	if strings.TrimSpace(string(out)) == "" {
		t.Error("expected at least a SAM header to be written")
	}

	mapContent, err := os.ReadFile(bcMap)
	if err != nil {
		t.Fatalf("ReadFile bcMap: %v", err)
	}
	if len(mapContent) != 0 {
		t.Errorf("expected empty bcMap for input with no records, got %q", mapContent)
	}
}
