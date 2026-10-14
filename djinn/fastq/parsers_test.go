package fastq

import (
	"os"
	"path/filepath"
	"testing"
)

func writeFastqFile(t *testing.T, lines ...string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "in.fq")
	var content string
	for _, l := range lines {
		content += l + "\n"
	}
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	return path
}

func TestCoreFqParser_Haplotagging(t *testing.T) {
	path := writeFastqFile(t,
		"@read1 BX:Z:A01C02B03D04", "ACGT", "+", "IIII",
	)
	parser, fmtType, err := CoreFqParser(path)
	if err != nil {
		t.Fatalf("CoreFqParser: %v", err)
	}
	if fmtType != "haplotagging" {
		t.Errorf("fmtType = %q, want %q", fmtType, "haplotagging")
	}
	if parser == nil {
		t.Fatal("expected a non-nil parser")
	}
}

func TestCoreFqParser_Stlfr(t *testing.T) {
	path := writeFastqFile(t,
		"@read1#123_456_789", "ACGT", "+", "IIII",
	)
	_, fmtType, err := CoreFqParser(path)
	if err != nil {
		t.Fatalf("CoreFqParser: %v", err)
	}
	if fmtType != "stlfr" {
		t.Errorf("fmtType = %q, want %q", fmtType, "stlfr")
	}
}

func TestCoreFqParser_Tellseq(t *testing.T) {
	path := writeFastqFile(t,
		"@read1:ATCGATCG", "ACGT", "+", "IIII",
	)
	_, fmtType, err := CoreFqParser(path)
	if err != nil {
		t.Fatalf("CoreFqParser: %v", err)
	}
	if fmtType != "tellseq" {
		t.Errorf("fmtType = %q, want %q", fmtType, "tellseq")
	}
}

func TestCoreFqParser_MixedFormatsIsError(t *testing.T) {
	path := writeFastqFile(t,
		"@read1 BX:Z:A01C02B03D04", "ACGT", "+", "IIII",
		"@read2#123_456_789", "ACGT", "+", "IIII",
	)
	_, _, err := CoreFqParser(path)
	if err == nil {
		t.Fatal("expected an error for mixed linked-read formats")
	}
}

func TestCoreFqParser_NoFormatFound(t *testing.T) {
	path := writeFastqFile(t,
		"@read1", "ACGT", "+", "IIII",
	)
	parser, fmtType, err := CoreFqParser(path)
	if err != nil {
		t.Fatalf("CoreFqParser: %v", err)
	}
	if parser != nil {
		t.Error("expected a nil parser when no format could be identified")
	}
	if fmtType != "" {
		t.Errorf("fmtType = %q, want empty", fmtType)
	}
}

func TestCoreFqParser_MissingFile(t *testing.T) {
	_, _, err := CoreFqParser("/nonexistent/path/does-not-exist.fq")
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}
