package xam

import (
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

// Returns a true if a barcode is valid (i.e. not invalid) in either
// haplotagging, tellseq, or stlfr formats
func IsValid(barcode string) bool {
	return !Invalid.MatchString(barcode)
}

// Convenience function to convert a boolean to integer,
// where false -> 0 and true -> 1.
func BoolToInt(vx bool) int {
	var i int
	if vx {
		i = 1
	}
	return i
}

// Convenience function to convert a boolean to integer string,
// where false -> "0" and true -> "1".
func BoolToSInt(vx bool) string {
	if vx {
		return "1"
	} else {
		return "0"
	}
}

// Search for an return a linked read barcode, whether a barcode
// was found (bool), and the value of the VX tag (bool). First searches for a BX tag,
// and if that isn't found, searches the record ID for a tellseq/stlfr style barcode.
// If nothing was found, returns ("", false, false). If a barcode was identified and
// a VX tag wasnt, the VX will be inferred from the barcode.
func FindBarcode(rec *sam.Record) (string, bool, bool) {
	bxVal, hasBX := GetStringTag(rec, "BX")
	vxVal, hasVX := GetVX(rec)
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
	if !hasVX && hasBX {
		vxVal = IsValid(bxVal)
	}
	return bxVal, hasBX, vxVal
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

/*
Convenience function to return the VX tag as a bool.
Returns (bool, bool), where the 2nd bool is whether
the VX tag was present. i.e.:

true, true -> VX:i = 1 and it was present

true, false -> VX:i = 1 and it wasn't present (impossible)

false, true -> VX:i = 0 and it was present

false, false VX:i = 0 because it wasn't present
*/
func GetVX(rec *sam.Record) (bool, bool) {
	vx, hasVX := GetIntTag(rec, "VX")
	vxVal := vx == 1
	return vxVal, hasVX
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

// Return the value of the string tag `tag` (XX:Z) for a sam.Record
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

// Return the value of the integer tag `tag` (XX:i) for a sam.Record
func GetIntTag(r *sam.Record, tag string) (int, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(int); ok {
				return s, true
			}
		}
	}
	return 0, false
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

// pairedFlag returns the SAM FLAG integer for a paired unmapped read.
// R1: 0x1|0x4|0x8|0x40 = 77   R2: 0x1|0x4|0x8|0x80 = 141
func PairedFlag(isRead1 bool) int {
	if isRead1 {
		return 77
	}
	return 141
}
