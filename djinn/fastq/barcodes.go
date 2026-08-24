package fastq

import (
	regexp "github.com/coregx/coregex"
)

var MissingBarcode []byte

var StlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)
var BxBarcode = regexp.MustCompile(`(?:BX\:Z\:(\S+))`)
var Invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")
var Tellseq = regexp.MustCompile(`:([ATCGN]+)(\s|$)`)

var Stlfr = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)(\s|$)`)
var StdBx = regexp.MustCompile(`BX:Z:(\S+)(?:\s|$)`)
var StdVx = regexp.MustCompile(`VX:i:([01])(?:\s|$)`)
var IlluminaOld = regexp.MustCompile(`/[12](?:\s|$)`)
var IlluminaNew = regexp.MustCompile(`[12]:[YN]:\d+:[A-Za-z0-9]+(?:\s|$)`)

var VXTAG = []byte{'V', 'X', ':', 'i', ':'}
var BXTAG = []byte{'B', 'X', ':', 'Z', ':'}
