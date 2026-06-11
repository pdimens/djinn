package xam

import (
	"log"
	"os"
	"regexp"
	"strconv"

	"github.com/biogo/hts/sam"
)

// Regex for invalid haplotagging, stlfr, tellseq barcodes
var Invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")

// Regex for valid haplotagging, stlfr, tellseq barcodes
var StlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)

// VX sam tag
var VxTag = sam.Tag{'V', 'X'}

// BX sam tag
var BxTag = sam.Tag{'B', 'X'}

// Returns a true if a barcode is invalid in either
// standard, haplotagging, tellseq, or stlfr formats
func IsInvalid(barcode string) bool {
	return Invalid.MatchString(barcode)
}

// Search for an return a linked read barcode and whether a barcode
// was found or not (bool). First searches for a BX tag, and if that isn't found, searches
// the record ID for a tellseq/stlfr style barcode. If nothing was found, returns
// ("", false)
func FindBarcode(rec *sam.Record) (string, bool) {
	bxVal, hasBX := GetStringTag(rec, "BX")
	if !hasBX {
		matches := StlfTell.FindStringSubmatch(rec.Name)
		if matches != nil {
			switch {
			case len(matches) > 1 && matches[1] != "":
				// matches[1] is the tellseq barcode e.g. "ATCGN"
				bxVal = matches[1]
				hasBX = true
			case len(matches) > 2 && matches[2] != "":
				// matches[2] is the stlfr barcode e.g. "1_2_3"
				bxVal = matches[2]
				hasBX = true
			}
		}
	}
	return bxVal, hasBX
}

// SetBX sets a string aux tag on a record
func SetBX(rec *sam.Record, val string) {
	aux := make(sam.Aux, 3+len(val)+1)
	aux[0] = 'B'
	aux[1] = 'X'
	aux[2] = 'Z'
	copy(aux[3:], val)
	aux[len(aux)-1] = 0 // null terminator

	for i, a := range rec.AuxFields {
		if a.Tag() == BxTag {
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
}

// SetVX sets an integer (0/1) auxiliary tag on a record.
func SetVX(rec *sam.Record, val bool) {
	b := byte(1)
	if val {
		b = 0
	}
	for i, a := range rec.AuxFields {
		if a.Tag() == VxTag {
			a[3] = b
			rec.AuxFields[i] = a
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, sam.Aux{'V', 'X', 'c', b})
}

// Return the value of the string tag `tag` for a sam.Record
func GetStringTag(r *sam.Record, tag string) (string, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(string); ok {
				return s, true
			}
		}
	}
	return "", false
}

// Create a copy of SAM header and add a new PG line for the djinn command `cmd`
func AddPG(src *sam.Header, cmd string) (*sam.Header, error) {
	// Clone the header via marshal/unmarshal
	b, err := src.MarshalText()
	if err != nil {
		return nil, err
	}
	dst := &sam.Header{}
	if err := dst.UnmarshalText(b); err != nil {
		return nil, err
	}

	// Build the @PG record
	pg := sam.NewProgram(
		"djinn",       // ID
		"djinn",       // name (PN)
		cmd,           // command line (CL)
		lastPGID(src), // previous PG ID (PP), or "" if none
		"1.0",         // version (VN) — set as appropriate
	)

	if err := dst.AddProgram(pg); err != nil {
		return nil, err
	}

	return dst, nil
}

// lastPGID returns the ID of the last @PG record in the header,
// which becomes the PP (previous program) of the new entry.
// Returns "" if there are no existing @PG records.
func lastPGID(h *sam.Header) string {
	progs := h.Progs()
	if len(progs) == 0 {
		return ""
	}
	return strconv.Itoa(progs[len(progs)-1].ID())
}

// Checks if the last argument is a file or input is coming from stdin
// returns either: file name, "-" (stdin), usage and exit.
func FileOrStdin(args []string, usage func()) string {
	var infile string
	switch len(args) {
	case 0:
		stat, err := os.Stdin.Stat()
		if err != nil {
			log.Fatalf("checking stdin: %v", err)
		}
		if (stat.Mode() & os.ModeCharDevice) != 0 {
			usage()
			os.Exit(1)
		}
		infile = "-"
	case 1:
		infile = args[0]
	default:
		usage()
		os.Exit(1)
	}
	return infile
}
