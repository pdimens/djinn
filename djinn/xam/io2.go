package xam

import (
	"bufio"
	"io"
	"log"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// Seqkit convenience function to check errors. If error != nil,
// does log.Fatalf(err)
func checkError(err error) {
	if err != nil {
		log.Fatalf("%v", err)
	}
}

// Create a new xam reader channel. Borrowed from seqkit. Usage:
//
// inChan, bamReader = NewBamReaderChan(inFile, chanCap, ioBuff, threads)
func NewXamReaderChan(inFile string, cp int, buff int, threads int) (chan *sam.Record, *AlignmentReader) {
	outChan := make(chan *sam.Record, cp)
	fh, err := os.Stdin, error(nil)
	if inFile != "-" {
		fh, err = os.Open(inFile)
		checkError(err)
	}

	bufReader := bufio.NewReaderSize(fh, buff)
	magic, err := bufReader.Peek(2)
	if err != nil {
		if fh != os.Stdin {
			fh.Close()
		}
		checkError(err)
	}

	var r AlignmentReader
	var sr *sam.Reader
	if magic[0] == 0x1f && magic[1] == 0x8b {
		r, err = bam.NewReader(bufReader, threads)
	} else {
		sr, err = sam.NewReader(bufReader)
		r = &SamReader{Reader: sr, closer: fh}
	}

	go func() {
		if fh != os.Stdin {
			defer fh.Close()
		}
		for {
			rec, err := r.Read()
			if err == io.EOF {
				close(outChan)
				return
			}
			if err != nil {
				close(outChan)
			}
			checkError(err)
			outChan <- rec
		}
	}()
	return outChan, &r
}

func NewXamSinkChan(cp int) (chan *sam.Record, chan bool) {
	outChan := make(chan *sam.Record, cp)
	doneChan := make(chan bool)
	go func() {
		for rec := range outChan {
			_ = rec
		}
		doneChan <- true
	}()

	return outChan, doneChan
}

// Create a sam/bam writer channel. Borrowed from seqkit. Usage:
//
// lastOut, doneChan = NewXamWriterChan(outFile, xamReader.Header(), chanCap, ioBuff, threads)
func NewXamWriterChan(outFile string, head *sam.Header, cp int, buff int, threads int, uncompressed bool) (chan *sam.Record, chan bool) {
	outChan := make(chan *sam.Record, cp)
	doneChan := make(chan bool)
	fh, err := os.Stdout, error(nil)
	if outFile != "-" {
		fh, err = os.Create(outFile)
		checkError(err)
	}

	var w AlignmentWriter
	bio := bufio.NewWriterSize(fh, buff)
	if uncompressed {
		sw, _ := sam.NewWriter(bio, head, sam.FlagDecimal)
		w = &SamWriter{Writer: sw}
	} else {
		w, err = bam.NewWriter(bio, head, threads)
	}

	go func() {
		for rec := range outChan {
			err := w.Write(rec)
			checkError(err)
		}
		w.Close()
		if !uncompressed {
			bio.Flush()
		}
		fh.Close()
		doneChan <- true
	}()
	return outChan, doneChan
}
