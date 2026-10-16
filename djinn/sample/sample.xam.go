package sample

import (
	"djinn/xam"
	"fmt"
	"math"
	"math/rand"
	"os"
	"slices"
	"time"
)

func SampleXam(infile string, downsample float64, seed, threads int, invalid, asSam bool) error {
	readThread, writeThread := xam.IoThreads(threads)

	// --- Get barcode map --------------------------
	infile = xam.FileOrStdin(infile)
	set := make(map[string]struct{}, 4_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── loop record channel to get barcodes ───────────────────────────────────
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" {
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
	// write alongside infile, not into the process's cwd: path.Base(infile)
	// dropped infile's directory entirely, writing wherever the command
	// happened to be run from.
	bcFile, err := os.Create(infile + ".bc")
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
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" {
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
