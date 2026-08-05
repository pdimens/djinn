package singletons

import (
	"djinn/xam"
	"io"
	"os"

	"github.com/biogo/hts/sam"
)

func getCount(infile string, threads int) map[string]int16 {
	infile = xam.FileOrStdin(infile)
	set := make(map[string]int16, 7_000_000)
	// ── open reader ───────────────────────────────────────────────────────────
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if hasBX && vxVal {
			set[bxVal]++
		}
	}
	return set
}

func FilterSingletons(infile, singletons string, asSam bool, threads int) error {
	err := xam.BamNotStdout(asSam)
	if err != nil {
		return err
	}
	// guard against draining stdin when getting barcode counts
	var f *os.File
	if infile == "-" {
		f, _ = os.CreateTemp("", "xam-spool-*.bam")
		io.Copy(f, os.Stdin)
		f.Seek(0, io.SeekStart)
		defer os.Remove(f.Name())
		infile = f.Name()
	}

	bcCounts := getCount(infile, threads)

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
	pg := xam.NewPG(hdr, "djinn sam filter-singletons "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)
	var singletonChan chan *sam.Record
	var singleDone chan bool

	if singletons != "" {
		singletonChan, singleDone = xam.NewXamWriterChan(singletons, hdr, xam.ChanCap, xam.IoBuf, 1, asSam)
	} else {
		singleDone = make(chan bool)
		close(singleDone)
	}

	// ── loop record channel ──────────────────────────────────────────────

	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if !hasBX || !vxVal {
			continue
		}
		if bcCounts[bxVal] > 2 {
			writeChan <- rec
		} else if singletonChan != nil {
			singletonChan <- rec
		}
	}
	close(writeChan)
	<-writeDone

	if singletonChan != nil {
		close(singletonChan)
		<-singleDone
	}
	return nil
}
