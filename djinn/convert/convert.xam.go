package convert

import (
	"djinn/barcodes"
	"djinn/xam"
	"fmt"
)

func ConvertXam(infile, convTo string, threads int, asSam bool) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator

	switch convTo {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
	case "stlfr":
		bcs = barcodes.NewStlfr()
	case "tellseq":
		bcs = barcodes.NewTellseq()
	case "10x":
		bcs = barcodes.NewTenX()
	default:
		return fmt.Errorf("unknown barcode type %q", convTo)
	}
	defer bcs.Close()

	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

	// ── open reader ───────────────────────────────────────────────────────────
	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam convert "+convTo+" "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel ──────────────────────────────────────────────
	var n int
	var ok bool
	seen := make(map[string][]byte, 4_000_000)

	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" || !vxVal {
			//TODO WRITE INVALID BARCODE
			bcs.InvalidInto(bcBuf)
			println("INVALID", string(bcBuf))
			continue
		}

		if converted, ok := seen[bxVal]; ok {
			// already seen — converted is the existing mapped value, no alloc happened for this lookup
			//TODO ADD OR REPLACE EXISTING BX:Z TAG WITH NEW BARCODE
			println(bxVal, string(converted))
			continue
		}

		// not seen — generate the next conversion barcode
		n, ok = bcs.NextInto(bcBuf)
		if !ok {
			return fmt.Errorf("Barcodes exhausted before conversion finished. The data provided has more unique barcodes than %q barcodes can support.", convTo)
		}
		converted := make([]byte, n) // must copy: buf gets reused next iteration
		copy(converted, bcBuf[:n])

		seen[bxVal] = converted
		println(bxVal, string(converted))
	}
	close(writeChan)
	<-writeDone

	return nil
}
