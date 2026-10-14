// Standalone fixture generator for benchmarking sort.SortByBX. Not part of
// the djinn build; run with `go run`. Writes an unsorted BAM file with
// random BX:Z tags, so the sort actually has work to do.
package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func main() {
	n := flag.Int("n", 3_500_000, "number of records")
	nBarcodes := flag.Int("barcodes", 200_000, "number of distinct BX values")
	seqLen := flag.Int("seqlen", 150, "read length")
	out := flag.String("out", "in.bam", "output BAM path")
	level := flag.Int("level", 0, "bam compression level (0 = none, matches sortAndSpill's temp files)")
	flag.Parse()

	hdr, err := sam.NewHeader(nil, []*sam.Reference{})
	if err != nil {
		log.Fatal(err)
	}
	ref, err := sam.NewReference("chr1", "", "", 1_000_000_000, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	if err := hdr.AddReference(ref); err != nil {
		log.Fatal(err)
	}

	f, err := os.Create(*out)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	bw, err := bam.NewWriterLevel(f, hdr, *level, 4)
	if err != nil {
		log.Fatal(err)
	}

	r := rand.New(rand.NewSource(1))
	bases := []byte("ACGT")
	seq := make([]byte, *seqLen)
	qual := make([]byte, *seqLen)
	for i := range qual {
		qual[i] = 30
	}
	barcodeBuf := make([]byte, 16)

	// Pre-generate the barcode pool once so repeats are exact (realistic:
	// many reads share a barcode) rather than every record getting a
	// fresh random string.
	barcodes := make([]string, *nBarcodes)
	for i := range barcodes {
		for j := range barcodeBuf {
			barcodeBuf[j] = bases[r.Intn(4)]
		}
		barcodes[i] = string(barcodeBuf)
	}

	log.Printf("generating %d records (%d distinct barcodes, %dbp, level=%d) -> %s", *n, *nBarcodes, *seqLen, *level, *out)
	for i := 0; i < *n; i++ {
		for j := range seq {
			seq[j] = bases[r.Intn(4)]
		}
		rec, err := sam.NewRecord(fmt.Sprintf("read%08d", i), nil, nil, -1, -1, 0, 0, nil, append([]byte(nil), seq...), append([]byte(nil), qual...), nil)
		if err != nil {
			log.Fatal(err)
		}
		rec.Flags = sam.Unmapped
		aux, err := sam.NewAux(sam.Tag{'B', 'X'}, barcodes[r.Intn(*nBarcodes)])
		if err != nil {
			log.Fatal(err)
		}
		rec.AuxFields = append(rec.AuxFields, aux)
		if err := bw.Write(rec); err != nil {
			log.Fatal(err)
		}
		if i%500_000 == 0 && i > 0 {
			log.Printf("  %d records written", i)
		}
	}
	if err := bw.Close(); err != nil {
		log.Fatal(err)
	}

	info, err := os.Stat(*out)
	if err == nil {
		log.Printf("done: %s is %.2f GB", *out, float64(info.Size())/1e9)
	}
}
