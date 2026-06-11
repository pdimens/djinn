package xam

import (
	"bufio"
	"io"
	"log"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── Constants ────────────────────────────────────────────────────────────────

const ChanCap = 2000
const IoBuf = 4 << 20 // 4 MiB

// ── Interfaces ────────────────────────────────────────────────────────────────

type AlignmentReader interface {
	Read() (*sam.Record, error)
	Header() *sam.Header
	Close() error
}

type AlignmentWriter interface {
	Write(*sam.Record) error
	Close() error
}

// ── SAM wrappers (sam.Reader/Writer lack Close()) ─────────────────────────────

type SamReader struct {
	*sam.Reader
	closer io.Closer
}

func (s *SamReader) Close() error {
	return s.closer.Close()
}

type SamWriter struct {
	*sam.Writer
}

func (s *SamWriter) Close() error {
	return nil
}

// ── Helpers ───────────────────────────────────────────────────────────────────

func checkError(err error) {
	if err != nil {
		log.Fatalf("%v", err)
	}
}

// FileOrStdin resolves a positional argument to a file path or "-" for stdin.
// Emits usage and exits if stdin is a terminal and no argument was provided.
func FileOrStdin(args []string, usage func()) string {
	switch len(args) {
	case 0:
		stat, err := os.Stdin.Stat()
		checkError(err)
		if (stat.Mode() & os.ModeCharDevice) != 0 {
			usage()
			os.Exit(1)
		}
		return "-"
	case 1:
		return args[0]
	default:
		usage()
		os.Exit(1)
	}
	return ""
}

// ── Reader channel ────────────────────────────────────────────────────────────

// NewXamReaderChan opens a SAM or BAM file (or stdin if inFile == "-") and
// streams records into the returned channel. The returned AlignmentReader
// exposes the file header. The goroutine owns all cleanup.
func NewXamReaderChan(inFile string, cp, buff, threads int) (chan *sam.Record, AlignmentReader) {
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

	// Build the reader. For SAM, fh ownership passes into SamReader so the
	// goroutine does not also close it. For BAM, the goroutine closes fh.
	var r AlignmentReader
	isBAM := magic[0] == 0x1f && magic[1] == 0x8b
	if isBAM {
		br, err := bam.NewReader(bufReader, threads)
		checkError(err)
		r = br
	} else {
		sr, err := sam.NewReader(bufReader)
		checkError(err)
		r = &SamReader{Reader: sr, closer: fh}
	}

	go func() {
		// Only close fh here for BAM; SamReader.Close() owns it for SAM.
		if isBAM && fh != os.Stdin {
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
				checkError(err) // fatal — exits, no fallthrough
				return          // defensive
			}
			outChan <- rec
		}
	}()

	return outChan, r
}

// ── Writer channel ────────────────────────────────────────────────────────────

// NewXamWriterChan writes records received on the returned channel to outFile
// (or stdout if outFile == "-"). Sends true on doneChan when the input channel
// is closed and all records have been flushed.
func NewXamWriterChan(outFile string, head *sam.Header, cp, buff, threads int, uncompressed bool) (chan *sam.Record, chan bool) {
	outChan := make(chan *sam.Record, cp)
	doneChan := make(chan bool)

	fh, err := os.Stdout, error(nil)
	if outFile != "-" {
		fh, err = os.Create(outFile)
		checkError(err)
	}

	bio := bufio.NewWriterSize(fh, buff)

	var w AlignmentWriter
	if uncompressed {
		sw, err := sam.NewWriter(bio, head, sam.FlagDecimal)
		checkError(err) // was silently discarded before
		w = &SamWriter{Writer: sw}
	} else {
		bw, err := bam.NewWriter(bio, head, threads)
		checkError(err)
		w = bw
	}

	go func() {
		for rec := range outChan {
			checkError(w.Write(rec))
		}
		w.Close()
		checkError(bio.Flush()) // always flush — was skipped for SAM before
		if fh != os.Stdout {
			fh.Close()
		}
		doneChan <- true
	}()

	return outChan, doneChan
}
