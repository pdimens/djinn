package xam

import (
	"djinn/barcodes"
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

func realTenXBarcodes(t *testing.T, n int) []string {
	t.Helper()
	g, err := barcodes.NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()
	buf := make([]byte, g.MaxLen())
	out := make([]string, 0, n)
	for i := 0; i < n; i++ {
		k, ok := g.NextInto(buf)
		if !ok {
			t.Fatalf("expected %d real 10X barcodes from the embedded whitelist", n)
		}
		out = append(out, string(buf[:k]))
	}
	return out
}

func TestDetectBarcodeFormat_Haplotagging(t *testing.T) {
	path := writeSamRecords(t, []string{"A01C02B03D04", "A05C06B07D08", "A09C10B11D12"})
	if got := DetectBarcodeFormat(path); got != "haplotagging" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "haplotagging")
	}
}

func TestDetectBarcodeFormat_Stlfr(t *testing.T) {
	path := writeSamRecords(t, []string{"1_2_3", "4_5_6", "7_8_9"})
	if got := DetectBarcodeFormat(path); got != "stlfr" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "stlfr")
	}
}

func TestDetectBarcodeFormat_Tellseq(t *testing.T) {
	// A nucleotide-shaped barcode that is not a real 10X whitelist entry.
	path := writeSamRecords(t, []string{"ATCGATCGATCGAAAA", "ATCGATCGATCGAAAC", "ATCGATCGATCGAAAG"})
	if got := DetectBarcodeFormat(path); got != "tellseq" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "tellseq")
	}
}

func TestDetectBarcodeFormat_TenX(t *testing.T) {
	path := writeSamRecords(t, realTenXBarcodes(t, 3))
	if got := DetectBarcodeFormat(path); got != "10x" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "10x")
	}
}

func TestDetectBarcodeFormat_MixedFormatsIsError(t *testing.T) {
	path := writeSamRecords(t, []string{"A01C02B03D04", "1_2_3"})
	if got := DetectBarcodeFormat(path); got != "error" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "error")
	}
}

func TestDetectBarcodeFormat_NoBarcodesIsError(t *testing.T) {
	path := writeSamRecords(t, []string{"", "", ""})
	if got := DetectBarcodeFormat(path); got != "error" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "error")
	}
}

func TestDetectBarcodeFormat_MissingFileIsError(t *testing.T) {
	if got := DetectBarcodeFormat("/nonexistent/path/does-not-exist.sam"); got != "error" {
		t.Errorf("DetectBarcodeFormat = %q, want %q", got, "error")
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
