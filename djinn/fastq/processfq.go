package fastq

import (
	"bytes"

	"github.com/shenwei356/bio/seqio/fastx"
)

const TabSep = '\t'
const Newline = '\n'
const PlusSign = '+'
const FastqAt = '@'

// Finds the BX and VX tags and removes the BX and VX tags
// along with the CASAVA /1 or 1:N:0:ATAG identifier. Returns
// a Barcode with the BX and VX tags. If VX isn't present, it's
// value is inferred from BX.
func HaplotagBX(rec *fastx.Record) bool {
	var valid bool

	// remove new CASAVA if present
	if iloc := IlluminaNew.FindIndex(rec.Desc); iloc != nil {
		rec.Desc = append(rec.Desc[:iloc[0]], rec.Desc[iloc[1]:]...)
	}

	// remove old CASAVA suffix if present
	if iloc := IlluminaOld.FindIndex(rec.ID); iloc != nil {
		rec.ID = append(rec.ID[:iloc[0]], rec.ID[iloc[1]:]...)
	}

	bloc := StdBx.FindSubmatchIndex(rec.Desc)
	if bloc == nil {
		if len(rec.Desc) > 0 {
			rec.Desc = append(rec.Desc, '\t')
		}
		rec.Desc = append(rec.Desc, MissingBarcode...)
		return false
	}
	bc := rec.Desc[bloc[2]:bloc[3]]
	vx := byte('0')

	vloc := StdVx.FindSubmatchIndex(rec.Desc)
	if vloc == nil {
		// no VX tag, check if barcode valid
		if !Invalid.Match(bc) {
			vx = '1'
			valid = true
		}
	} else if bytes.Equal(rec.Desc[vloc[2]:vloc[3]], []byte{'1'}) {
		// VX tag present and is 1
		vx = '1'
		valid = true
	}

	// remove BX/VX from rec.Desc: later match first, so the earlier
	// match's indices stay valid and we never shift bytes we're
	// about to discard
	if vloc != nil {
		first, second := bloc, vloc
		if vloc[0] < bloc[0] {
			first, second = vloc, bloc
		}
		rec.Desc = append(rec.Desc[:second[0]], rec.Desc[second[1]:]...)
		rec.Desc = append(rec.Desc[:first[0]], rec.Desc[first[1]:]...)
	} else {
		rec.Desc = append(rec.Desc[:bloc[0]], rec.Desc[bloc[1]:]...)
	}

	// add BX and VX back in
	if len(rec.Desc) > 0 {
		rec.Desc = append(rec.Desc, '\t')
	}
	rec.Desc = append(rec.Desc, VXTAG...)
	rec.Desc = append(rec.Desc, vx, '\t')
	rec.Desc = append(rec.Desc, BXTAG...)
	rec.Desc = append(rec.Desc, bc...)

	return valid
}

// Find tellseq barcode inline in rec.ID, checking it for validity.
// Converts rec in-place into standard format. Returns early if no barcode was found.
// Returns a bool of whether the barcode was valid (true) or not (false) as a sentinel
// value for how to post-process the read.
func Tellseq2Std(rec *fastx.Record) bool {
	var valid bool

	// remove new CASAVA if present
	if iloc := IlluminaNew.FindIndex(rec.Desc); iloc != nil {
		rec.Desc = append(rec.Desc[:iloc[0]], rec.Desc[iloc[1]:]...)
	}

	// remove old CASAVA suffix if present
	if iloc := IlluminaOld.FindIndex(rec.ID); iloc != nil {
		rec.ID = append(rec.ID[:iloc[0]], rec.ID[iloc[1]:]...)
	}

	// find stlfr barcode in the record ID
	bloc := Tellseq.FindSubmatchIndex(rec.ID)
	if bloc == nil {
		if len(rec.Desc) > 0 {
			rec.Desc = append(rec.Desc, '\t')
		}
		rec.Desc = append(rec.Desc, MissingBarcode...)
		return false
	}

	// bc aliases rec.ID's backing array — must be read before the
	// in-place shift below overwrites this region
	bc := rec.ID[bloc[2]:bloc[3]]
	vx := byte('0')
	if !Invalid.Match(bc) {
		valid = true
		vx = '1'
	}

	// write tag directly into rec.Desc — no intermediate buffer
	if len(rec.Desc) > 0 {
		rec.Desc = append(rec.Desc, '\t')
	}
	rec.Desc = append(rec.Desc, VXTAG...)
	rec.Desc = append(rec.Desc, vx, '\t')
	rec.Desc = append(rec.Desc, BXTAG...)
	rec.Desc = append(rec.Desc, bc...)

	// remove the full match from rec.ID in place
	rec.ID = append(rec.ID[:bloc[0]], rec.ID[bloc[1]:]...)

	return valid
}

// Find tellseq barcode inline in rec.ID, checking it for validity.
// Converts rec in-place into standard format. Returns early if no barcode was found.
// Returns a bool of whether the barcode was valid (true) or not (false) as a sentinel
// value for how to post-process the read.
func Stlfr2Std(rec *fastx.Record) bool {
	var valid bool

	// remove new CASAVA if present
	if iloc := IlluminaNew.FindIndex(rec.Desc); iloc != nil {
		rec.Desc = append(rec.Desc[:iloc[0]], rec.Desc[iloc[1]:]...)
	}

	// remove old CASAVA suffix if present
	if iloc := IlluminaOld.FindIndex(rec.ID); iloc != nil {
		rec.ID = append(rec.ID[:iloc[0]], rec.ID[iloc[1]:]...)
	}

	// find stlfr barcode
	bloc := Stlfr.FindSubmatchIndex(rec.ID)
	if bloc == nil {
		if len(rec.Desc) > 0 {
			rec.Desc = append(rec.Desc, '\t')
		}
		rec.Desc = append(rec.Desc, MissingBarcode...)
		return false
	}

	// bc aliases rec.ID's backing array — must be read before the
	// in-place shift below overwrites this region
	bc := rec.ID[bloc[2]:bloc[3]]
	vx := byte('0')
	if !Invalid.Match(bc) {
		valid = true
		vx = '1'
	}

	// write tag directly into rec.Desc — no intermediate buffer
	if len(rec.Desc) > 0 {
		rec.Desc = append(rec.Desc, '\t')
	}
	rec.Desc = append(rec.Desc, VXTAG...)
	rec.Desc = append(rec.Desc, vx, '\t')
	rec.Desc = append(rec.Desc, BXTAG...)
	rec.Desc = append(rec.Desc, bc...)

	// remove the full match from rec.ID in place
	rec.ID = append(rec.ID[:bloc[0]], rec.ID[bloc[1]:]...)

	return valid
}
