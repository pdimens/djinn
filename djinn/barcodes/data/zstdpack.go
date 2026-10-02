//go:build ignore

// Standalone helper, not part of the djinn build (see the ignore tag
// above): recompresses a barcode-list text file -- one barcode per line,
// optionally itself gzip'd -- into zstd at maximum effort, for embedding
// via go:embed as 10x_whitelist.txt.zst. See README.md in this directory.
//
// Run it with `go run` from inside the djinn/djinn module (so it resolves
// github.com/klauspost/compress from the existing go.sum):
//
//	go run barcodes/data/zstdpack.go -in whitelist.txt.gz -out barcodes/data/10x_whitelist.txt.zst
package main

import (
	"compress/gzip"
	"flag"
	"io"
	"log"
	"os"
	"strings"

	"github.com/klauspost/compress/zstd"
)

func main() {
	in := flag.String("in", "", "input barcode list (.txt or .txt.gz)")
	out := flag.String("out", "", "output .zst path")
	flag.Parse()
	if *in == "" || *out == "" {
		log.Fatal("usage: zstdpack -in <whitelist.txt[.gz]> -out <whitelist.txt.zst>")
	}

	f, err := os.Open(*in)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(*in, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			log.Fatal(err)
		}
		defer gz.Close()
		r = gz
	}

	outF, err := os.Create(*out)
	if err != nil {
		log.Fatal(err)
	}
	defer outF.Close()

	enc, err := zstd.NewWriter(outF, zstd.WithEncoderLevel(zstd.SpeedBestCompression))
	if err != nil {
		log.Fatal(err)
	}
	n, err := io.Copy(enc, r)
	if err != nil {
		log.Fatal(err)
	}
	if err := enc.Close(); err != nil {
		log.Fatal(err)
	}

	outInfo, statErr := os.Stat(*out)
	if statErr != nil {
		log.Fatal(statErr)
	}
	log.Printf("wrote %s: %d bytes in -> %d bytes out (%.2fx)", *out, n, outInfo.Size(), float64(n)/float64(outInfo.Size()))
}
