package xam

import (
	"encoding/binary"
	"strings"

	"github.com/biogo/hts/sam"
	regexp "github.com/coregx/coregex"
)

// Regex for invalid haplotagging, stlfr, tellseq barcodes
//var Invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")

// Regex for valid haplotagging, stlfr, tellseq barcodes
var StlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)

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

//func IsValid(barcode string) bool {
//	return !Invalid.MatchString(barcode)
//}

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

// Search for an return a linked read barcode, ("" if not found) and the value of the VX tag (bool).
// First searches for a BX tag, and if that isn't found, searches the record ID for a tellseq/stlfr style barcode.
// If nothing was found, returns ("", false). If a barcode was identified and
// a VX tag wasnt, the VX will be inferred from the barcode.
func FindBarcode(rec *sam.Record) (string, bool) {
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
	if !hasVX && bxVal != "" {
		vxVal = IsValid(bxVal)
	}
	return bxVal, vxVal
}

// SetBX sets a string aux tag on a record
func SetBxString(rec *sam.Record, val *string) {
	//aux := make(sam.Aux, 3+len(val)+1)
	aux := make(sam.Aux, 3+len(*val))
	aux[0], aux[1], aux[2] = 'B', 'X', 'Z'
	copy(aux[3:], *val)
	// aux[len(aux)-1] = 0 // null terminator

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
	//aux := make(sam.Aux, 3+len(val)+1)
	aux := make(sam.Aux, 3+len(*val))
	aux[0], aux[1], aux[2] = 'B', 'X', 'Z'
	copy(aux[3:], *val)
	// aux[len(aux)-1] = 0 // null terminator

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

	for _, a := range rec.AuxFields {
		if a.Tag() == VxTag {
			copy(a[3:7], buf)
			return
		}
	}
	aux := sam.Aux{'V', 'X', 'i', 0, 0, 0, 0}
	copy(aux[3:7], buf)
	rec.AuxFields = append(rec.AuxFields, aux)
}

/*

func SetVX(rec *sam.Record, isValid bool) {
	b := byte(0)
	if isValid {
		b = 1
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
*/

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

// pairedFlag returns the SAM FLAG integer for a paired unmapped read.
// R1: 0x1|0x4|0x8|0x40 = 77   R2: 0x1|0x4|0x8|0x80 = 141
func PairedFlag(isRead1 bool) int {
	if isRead1 {
		return 77
	}
	return 141
}
