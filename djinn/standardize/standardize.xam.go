package standardize

import "djinn/xam"

func Standardize(infile string, threads int, asSam bool) error {
	infile = xam.FileOrStdin(infile)
	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, r := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := r.Header()
	pg := xam.NewPG(hdr, "djinn standardize "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	outChan, doneChan := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel ──────────────────────────────────────────────

	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if hasBX {
			xam.SetBX(rec, bxVal)
			xam.SetVX(rec, vxVal)
		}
		//fmt.Println(bxVal, vxVal)
		// push updated record into writer channel
		outChan <- rec
	}

	close(outChan)

	<-doneChan
	return nil
}
