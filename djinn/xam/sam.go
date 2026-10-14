package xam

import (
	"encoding/binary"
	"strings"

	"github.com/biogo/hts/sam"
)

// VX sam tag
var VxTag = sam.Tag{'V', 'X'}

// BX sam tag
var BxTag = sam.Tag{'B', 'X'}

// Returns a true if a barcode is valid (i.e. not invalid) in either
// haplotagging, tellseq, or stlfr formats
func IsValid(barcode string) bool {
	if strings.IndexByte(barcode, 'N') != -1 {
		return false
	}
	if strings.HasPrefix(barcode, "0_") || strings.HasSuffix(barcode, "_0") || strings.Contains(barcode, "_0_") {
		return false
	}
	for _, p := range [4]string{"A00", "B00", "C00", "D00"} {
		if strings.Contains(barcode, p) {
			return false
		}
	}
	return true
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

// isBase reports whether b is one of the DNA bases used in barcode text
// (including the ambiguity code N).
func isBase(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'N':
		return true
	}
	return false
}

func isDigit(b byte) bool { return b >= '0' && b <= '9' }

// isStlfrBody reports whether s matches "[0-9]+_[0-9]+_[0-9]+" exactly.
func isStlfrBody(s string) bool {
	parts, start, digits := 0, 0, false
	for i := 0; i <= len(s); i++ {
		if i == len(s) || s[i] == '_' {
			if !digits || i == start {
				return false
			}
			parts++
			start = i + 1
			digits = false
			continue
		}
		if !isDigit(s[i]) {
			return false
		}
		digits = true
	}
	return parts == 3
}

// findInlineBarcode locates a linked-read barcode inline in name, in
// either TELLseq (":<bases>" suffix) or stLFR ("#<n>_<n>_<n>" suffix)
// style, mirroring what the pattern `(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`
// matched, but as a pair of native byte scans instead of a regexp.
func findInlineBarcode(name string) (barcode string, ok bool) {
	// TELLseq: a run of bases at the very end of name, preceded by ':'.
	n := len(name)
	i := n
	for i > 0 && isBase(name[i-1]) {
		i--
	}
	if i < n && i > 0 && name[i-1] == ':' {
		return name[i:n], true
	}

	// stLFR: "#<n>_<n>_<n>" at the end of name. Scan every '#' left to
	// right (mirroring the leftmost-match behavior of the regexp) and
	// take the first whose remainder matches the full n_n_n shape.
	for j := 0; j < n; j++ {
		if name[j] != '#' {
			continue
		}
		if body := name[j+1:]; isStlfrBody(body) {
			return body, true
		}
	}
	return "", false
}

// Search for an return a linked read barcode, ("" if not found) and the value of the VX tag (bool).
// First searches for a BX tag, and if that isn't found, searches the record ID for a tellseq/stlfr style barcode.
// If nothing was found, returns ("", false). If a barcode was identified and
// a VX tag wasnt, the VX will be inferred from the barcode.
func FindBarcode(rec *sam.Record) (string, bool) {
	bxVal, hasBX := GetStringTag(rec, "BX")
	vxVal, hasVX := GetVX(rec)
	if !hasBX {
		bxVal, hasBX = findInlineBarcode(rec.Name)
	}
	if !hasVX && bxVal != "" {
		vxVal = IsValid(bxVal)
	}
	return bxVal, vxVal
}

// SetBX sets a string aux tag on a record
func SetBxString(rec *sam.Record, val *string) {
	aux := make(sam.Aux, 3+len(*val))
	aux[0], aux[1], aux[2] = 'B', 'X', 'Z'
	copy(aux[3:], *val)

	for i, a := range rec.AuxFields {
		if a.Tag() == BxTag {
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
}

// SetBX sets a string aux tag on a record
func SetBxByte(rec *sam.Record, val *[]byte) {
	aux := make(sam.Aux, 3+len(*val))
	aux[0], aux[1], aux[2] = 'B', 'X', 'Z'
	copy(aux[3:], *val)

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
func SetVX(rec *sam.Record, isValid bool) {
	var v int32
	if isValid {
		v = 1
	}
	buf := make([]byte, 4)
	binary.LittleEndian.PutUint32(buf, uint32(v))

	aux := sam.Aux{'V', 'X', 'i', 0, 0, 0, 0}
	copy(aux[3:7], buf)

	for i, a := range rec.AuxFields {
		if a.Tag() == VxTag {
			// Replace the whole field rather than copying into the
			// existing byte range: an aux field parsed from SAM text
			// (e.g. "VX:i:1") may have been encoded by the library as
			// a narrower type ('c'/'C'/'s'/'S', 1-2 bytes) instead of
			// the 4-byte 'i' this function assumes, and copy(a[3:7], buf)
			// would panic with an out-of-range slice on such a field.
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
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
			// sam.Aux.Value() never returns a plain "int" for integer-kind
			// tags: depending on the field's stored width it returns
			// int8/uint8/int16/uint16/int32/uint32, so a bare .(int)
			// assertion never succeeds and this always reported "not found".
			switch v := aux.Value().(type) {
			case int8:
				return int(v), true
			case uint8:
				return int(v), true
			case int16:
				return int(v), true
			case uint16:
				return int(v), true
			case int32:
				return int(v), true
			case uint32:
				return int(v), true
			}
		}
	}
	return 0, false
}

// pairedFlag returns the SAM FLAG integer for a paired unmapped read.
// R1: 0x1|0x4|0x8|0x40 = 77   R2: 0x1|0x4|0x8|0x80 = 141
func PairedFlag(isRead1 bool) int {
	if isRead1 {
		return 77
	}
	return 141
}
