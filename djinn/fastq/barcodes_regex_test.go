package fastq

import "testing"

func TestInvalidRegex(t *testing.T) {
	tests := []struct {
		name string
		in   string
		want bool
	}{
		{"contains N", "ATCGN", true},
		{"haplotag A00", "A00C12B34D56", true},
		{"haplotag B00", "A12C12B00D56", true},
		{"haplotag D00 at end", "A12C12B34D00", true},
		{"stlfr leading zero", "0_5_5", true},
		{"stlfr middle zero", "5_0_5", true},
		{"stlfr trailing zero", "5_5_0", true},
		{"all valid haplotag", "A12C34B56D78", false},
		{"all valid stlfr", "5_5_5", false},
		{"pure ATCG", "ATCG", false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := Invalid.Match([]byte(tt.in))
			if got != tt.want {
				t.Errorf("Invalid.Match(%q) = %v, want %v", tt.in, got, tt.want)
			}
		})
	}
}

func TestStlfrRegex(t *testing.T) {
	tests := []struct {
		in      string
		wantHit bool
		want    string
	}{
		{"read1#123_456_789", true, "123_456_789"},
		{"read1#123_456_789 comment", true, "123_456_789"},
		{"read1", false, ""},
		{"read1#123_456", false, ""},
	}
	for _, tt := range tests {
		loc := Stlfr.FindSubmatchIndex([]byte(tt.in))
		if (loc != nil) != tt.wantHit {
			t.Errorf("Stlfr.FindSubmatchIndex(%q) hit=%v, want %v", tt.in, loc != nil, tt.wantHit)
			continue
		}
		if tt.wantHit {
			got := tt.in[loc[2]:loc[3]]
			if got != tt.want {
				t.Errorf("Stlfr match = %q, want %q", got, tt.want)
			}
		}
	}
}

func TestTellseqRegex(t *testing.T) {
	tests := []struct {
		in      string
		wantHit bool
		want    string
	}{
		{"read1:ATCGATCG", true, "ATCGATCG"},
		{"read1:ATCGATCG comment", true, "ATCGATCG"},
		{"read1", false, ""},
	}
	for _, tt := range tests {
		loc := Tellseq.FindSubmatchIndex([]byte(tt.in))
		if (loc != nil) != tt.wantHit {
			t.Errorf("Tellseq.FindSubmatchIndex(%q) hit=%v, want %v", tt.in, loc != nil, tt.wantHit)
			continue
		}
		if tt.wantHit {
			got := tt.in[loc[2]:loc[3]]
			if got != tt.want {
				t.Errorf("Tellseq match = %q, want %q", got, tt.want)
			}
		}
	}
}

func TestStdBxRegex(t *testing.T) {
	tests := []struct {
		in      string
		wantHit bool
		want    string
	}{
		{"BX:Z:A01C02B03D04", true, "A01C02B03D04"},
		{"VX:i:1\tBX:Z:A01C02B03D04", true, "A01C02B03D04"},
		{"no tag here", false, ""},
	}
	for _, tt := range tests {
		loc := StdBx.FindSubmatchIndex([]byte(tt.in))
		if (loc != nil) != tt.wantHit {
			t.Errorf("StdBx.FindSubmatchIndex(%q) hit=%v, want %v", tt.in, loc != nil, tt.wantHit)
			continue
		}
		if tt.wantHit {
			got := tt.in[loc[2]:loc[3]]
			if got != tt.want {
				t.Errorf("StdBx match = %q, want %q", got, tt.want)
			}
		}
	}
}

func TestStdVxRegex(t *testing.T) {
	tests := []struct {
		in      string
		wantHit bool
		want    string
	}{
		{"VX:i:1", true, "1"},
		{"VX:i:0", true, "0"},
		{"VX:i:2", false, ""},
		{"no tag", false, ""},
	}
	for _, tt := range tests {
		loc := StdVx.FindSubmatchIndex([]byte(tt.in))
		if (loc != nil) != tt.wantHit {
			t.Errorf("StdVx.FindSubmatchIndex(%q) hit=%v, want %v", tt.in, loc != nil, tt.wantHit)
			continue
		}
		if tt.wantHit {
			got := tt.in[loc[2]:loc[3]]
			if got != tt.want {
				t.Errorf("StdVx match = %q, want %q", got, tt.want)
			}
		}
	}
}

func TestIlluminaOldRegex(t *testing.T) {
	tests := []struct {
		in      string
		wantHit bool
	}{
		{"read1/1", true},
		{"read1/2", true},
		{"read1/3", false},
		{"read1", false},
	}
	for _, tt := range tests {
		got := IlluminaOld.Match([]byte(tt.in)) != false && IlluminaOld.FindIndex([]byte(tt.in)) != nil
		if got != tt.wantHit {
			t.Errorf("IlluminaOld match(%q) = %v, want %v", tt.in, got, tt.wantHit)
		}
	}
}

func TestIlluminaNewRegex(t *testing.T) {
	tests := []struct {
		in      string
		wantHit bool
	}{
		{"1:N:0:ATCGATCG", true},
		{"2:Y:0:ATCGATCG", true},
		{"3:N:0:ATCGATCG", false},
		{"no casava here", false},
	}
	for _, tt := range tests {
		got := IlluminaNew.FindIndex([]byte(tt.in)) != nil
		if got != tt.wantHit {
			t.Errorf("IlluminaNew match(%q) = %v, want %v", tt.in, got, tt.wantHit)
		}
	}
}
