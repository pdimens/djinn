package concat

import (
	"djinn/barcodes"
	"djinn/xam"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func ConcatXam(infiles []string, threads int, asSam bool) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator
	var err error

	bcFmt, err := xam.DetectBarcodeFormat(infiles[0])
	if err != nil {
		return err
	}

	switch bcFmt {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
	case "stlfr":
		bcs = barcodes.NewStlfr()
	case "nucleotide-18":
		bcs = barcodes.NewTellseq()
	case "nucleotide-16":
		bcs = barcodes.NewTenX()
	default:
		// generic-XX
		n, err := strconv.Atoi(bcFmt[len(bcFmt)-2:])
		if err != nil {
			return err
		}
		bcs, err = barcodes.NewGeneric(n)
		if err != nil {
			return err
		}
	}
	defer bcs.Close()

	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}

	// ── open first file for reading────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infiles[0], xam.ChanCap, xam.IoBuf, readThread)

	// update PG line in header
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam concat "+strings.Join(os.Args[1:], " "))
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// scale to number of input files
	seenGlobal := make(map[string]struct{}, 1_000_000*len(infiles))

	// ── loop records of first file ──────────────────────────────────────────────
	// this will populate the globally-seen map
	invalidBarcode := bcs.GetInvalid()
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" || !vxVal {
			xam.SetBxByte(rec, &invalidBarcode)
			writeChan <- rec
			continue
		}
		if _, ok := seenGlobal[bxVal]; !ok {
			seenGlobal[bxVal] = struct{}{}
		}
		writeChan <- rec
	}
	br.Close()

	// ── loop through remaining files ──────────────────────────────────────────────
	var n int
	var _ok bool

	for _, samfile := range infiles[1:] {
		recChan, br := xam.NewXamReaderChan(samfile, xam.ChanCap, xam.IoBuf, readThread)
		seenLocal := make(map[string][]byte, 2_000_000)

		for rec := range recChan {
			bxVal, vxVal := xam.FindBarcode(rec)
			if bxVal == "" || !vxVal {
				xam.SetBxByte(rec, &invalidBarcode)
				writeChan <- rec
				continue
			}
			// is it seen in the global inventory?
			_, ok := seenGlobal[bxVal]
			switch ok {
			case false:
				// not previously seen in any file, add to global inventory
				seenGlobal[bxVal] = struct{}{}
				// can keep this barcode, no need for conversion, add to local as self
				seenLocal[bxVal] = []byte(bxVal)
			case true:
				// previously seen in another file, add to local inventory with new barcode
				n, _ok = bcs.NextInto(bcBuf)
				if !_ok {
					return fmt.Errorf("The data provided has more unique barcodes than the input technology (%q) can support. Stopped at file %q.", bcFmt, samfile)
				}
				converted := make([]byte, n) // must copy: buf gets reused next iteration
				copy(converted, bcBuf[:n])
				seenLocal[bxVal] = converted

				xam.SetBxByte(rec, &converted)
				writeChan <- rec
			}
			writeChan <- rec
		}
		br.Close()
	}
	close(writeChan)
	<-writeDone
	return nil
}
