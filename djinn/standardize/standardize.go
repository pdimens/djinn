package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"strings"

	"djinn/xam"
)

// TODO need a mechanism for barcode conversion
func main() {
	nThreads := flag.Int("threads", 1, "Number of threads for BAM io. Diminishing returns beyond 4.")
	asSam := flag.Bool("sam", false, "Output as uncompressed SAM")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: djinn sam standardize [options] [bam|sam]\n\n")
		fmt.Fprintf(os.Stderr, "options:\n")
		flag.PrintDefaults()
	}
	flag.Parse()

	args := flag.Args()
	infile := xam.FileOrStdin(args, flag.Usage)

	threads := min(max(*nThreads, 1), runtime.NumCPU())
	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, r := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)
	newHeader, err := xam.AddPG(r.Header(), strings.Join(os.Args, " "))
	if err != nil {
		log.Fatal(err)
	}

	if !*asSam {
		stat, _ := os.Stdout.Stat()
		if (stat.Mode() & os.ModeCharDevice) != 0 {
			log.Fatal("refusing to write BAM to the terminal, use -sam or pipe to another tool/file")
		}
	}
	// ── open writer ───────────────────────────────────────────────────────────
	outChan, doneChan := xam.NewXamWriterChan(
		"-", newHeader, xam.ChanCap, xam.IoBuf, writeThread, *asSam,
	)

	// ── loop record channel ──────────────────────────────────────────────

	for rec := range recChan {
		bxVal, hasBX, VX := xam.FindBarcode(rec)
		if hasBX {
			xam.SetBX(rec, bxVal)
			xam.SetVX(rec, VX)
		}
		// push updated record into writer channel
		outChan <- rec
	}

	close(outChan) // signal writer goroutine that we're done

	<-doneChan // wait for writer to flush and close
}
