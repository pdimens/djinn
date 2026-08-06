package sample

import (
	"djinn/xam"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path"
	"slices"
	"time"
)

func Sample(infile string, downsample float64, seed, threads int, invalid, asSam bool) error {
	err := xam.BamNotStdout(asSam)
	if err != nil {
		return err
	}
	// --- Get barcode map --------------------------
	infile = xam.FileOrStdin(infile)
	set := make(map[string]struct{}, 7_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── loop record channel to get barcodes ───────────────────────────────────
	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if !hasBX {
			continue
		}
		if !vxVal && !invalid {
			continue
		}
		if _, ok := set[bxVal]; ok {
			continue
		} else {
			set[bxVal] = struct{}{}
		}
	}

	// ---- Sort and extract barcodes to keep -----------------------------------
	nBC := len(set)

	if float64(nBC) < downsample {
		return fmt.Errorf("The input has fewer barcodes (%v) than the requested downsampling amount (%v)", nBC, downsample)
	}

	keys := make([]string, 0, nBC)
	for k := range set {
		keys = append(keys, k)
	}
	slices.Sort(keys) // ensures reproducibility later

	// downsample by shuffling barcodes -> selecting first n
	var r *rand.Rand
	if seed >= 0 {
		r = rand.New(rand.NewSource(int64(seed)))
	} else {
		r = rand.New(rand.NewSource(time.Now().UTC().UnixNano()))
	}
	r.Shuffle(nBC, func(i, j int) { keys[i], keys[j] = keys[j], keys[i] })

	var nBcToKeep = int(downsample)
	bcFile, err := os.Create(path.Base(infile) + ".bc")
	if err != nil {
		return err
	}

	if downsample < 1.0 {
		nBcToKeep = int(math.Round(float64(nBC) * downsample))
	}
	bcToKeep := make(map[string]struct{}, nBcToKeep)
	for i := range nBcToKeep {
		bcToKeep[keys[i]] = struct{}{}
		bcFile.WriteString(keys[i] + "\n")
	}

	err = bcFile.Close()
	if err != nil {
		return err
	}

	// ── open reader (again) ───────────────────────────────────────────────────
	recChan, br = xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam sample "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}
	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel to get barcodes ───────────────────────────────────
	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if !hasBX {
			continue
		}
		if !vxVal && !invalid {
			continue
		}
		if _, ok := bcToKeep[bxVal]; ok {
			writeChan <- rec
		}
	}
	close(writeChan)
	<-writeDone

	return nil
}
