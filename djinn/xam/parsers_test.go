package xam

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const parsersTestHeader = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"

func writeSamRecords(t *testing.T, bxValues []string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "in.sam")

	var sb strings.Builder
	sb.WriteString(parsersTestHeader)
	for i, bx := range bxValues {
		sb.WriteString(fmt.Sprintf("read%d\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", i))
		if bx != "" {
			sb.WriteString("\tBX:Z:")
			sb.WriteString(bx)
		}
		sb.WriteByte('\n')
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	return path
}

func TestDetectBarcodeFormat_Haplotagging(t *testing.T) {
	path := writeSamRecords(t, []string{"A01C02B03D04", "A05C06B07D08", "A09C10B11D12"})
	if got, _ := DetectBarcodeFormat(path); got != "haplotagging" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "haplotagging")
	}
}

func TestDetectBarcodeFormat_Stlfr(t *testing.T) {
	path := writeSamRecords(t, []string{"1_2_3", "4_5_6", "7_8_9"})
	if got, _ := DetectBarcodeFormat(path); got != "stlfr" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "stlfr")
	}
}

func TestDetectBarcodeFormat_GenericTellseqLength(t *testing.T) {
	// 18bp nucleotide barcodes -- TELLseq length, reported generically.
	path := writeSamRecords(t, []string{"ATCGATCGATCGATCGAA", "ATCGATCGATCGATCGAC", "ATCGATCGATCGATCGAG"})
	if got, _ := DetectBarcodeFormat(path); got != "generic-18" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "generic-18")
	}
}

func TestDetectBarcodeFormat_GenericTenXLength(t *testing.T) {
	// 16bp nucleotide barcodes -- 10X length, reported generically.
	path := writeSamRecords(t, []string{"ATCGATCGATCGAAAA", "ATCGATCGATCGAAAC", "ATCGATCGATCGAAAG"})
	if got, _ := DetectBarcodeFormat(path); got != "generic-16" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "generic-16")
	}
}

func TestDetectBarcodeFormat_MixedFormatsIsError(t *testing.T) {
	path := writeSamRecords(t, []string{"A01C02B03D04", "1_2_3"})
	if _, err := DetectBarcodeFormat(path); err == nil {
		t.Errorf("DetectBarcodeFormat = %q, want %q", err, "error")
	}
}

func TestDetectBarcodeFormat_NoBarcodesIsError(t *testing.T) {
	path := writeSamRecords(t, []string{"", "", ""})
	if _, err := DetectBarcodeFormat(path); err == nil {
		t.Errorf("DetectBarcodeFormat = %q, want %q", err, "error")
	}
}

func TestDetectBarcodeFormat_MissingFileIsError(t *testing.T) {
	if _, err := DetectBarcodeFormat("/nonexistent/path/does-not-exist.sam"); err == nil {
		t.Errorf("DetectBarcodeFormat = %q, want %q", err, "error")
	}
}

func TestIsHaplotaggingBarcode(t *testing.T) {
	cases := []struct {
		bc   string
		want bool
	}{
		{"A01C02B03D04", true},
		{"A99C99B99D99", true},
		{"a01C02B03D04", false}, // lowercase leading letter
		{"A01C02B03D0", false},  // too short
		{"A01C02B03D041", false},
		{"ATCGATCGATCG", false}, // nucleotide, wrong shape
		{"1_2_3", false},
	}
	for _, c := range cases {
		if got := isHaplotaggingBarcode(c.bc); got != c.want {
			t.Errorf("isHaplotaggingBarcode(%q) = %v, want %v", c.bc, got, c.want)
		}
	}
}

func TestIsNucleotideBarcode(t *testing.T) {
	cases := []struct {
		bc   string
		want bool
	}{
		{"ATCGN", true},
		{"ATCG", true},
		{"", false},
		{"A01C02B03D04", false},
		{"1_2_3", false},
		{"ATCGX", false},
	}
	for _, c := range cases {
		if got := isNucleotideBarcode(c.bc); got != c.want {
			t.Errorf("isNucleotideBarcode(%q) = %v, want %v", c.bc, got, c.want)
		}
	}
}
