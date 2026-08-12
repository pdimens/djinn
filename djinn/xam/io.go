package xam

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── Constants ────────────────────────────────────────────────────────────────

const ChanCap = 5000
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

func BamNotStdout(isSam bool) error {
	stat, err := os.Stdout.Stat()
	if err != nil {
		return err
	}
	if (stat.Mode()&os.ModeCharDevice) != 0 && !isSam {
		return fmt.Errorf("refusing to write BAM to the terminal, as this was likely unintended. Use -S or pipe to another tool/file")
	}
	return nil
}

// FileOrStdin resolves a positional argument to a file path or "-" for stdin.
func FileOrStdin(infile string) string {
	if len(infile) == 0 {
		stat, err := os.Stdin.Stat()
		checkError(err)
		if (stat.Mode() & os.ModeCharDevice) != 0 {
			os.Exit(1)
		}
		return "-"
	}
	return infile
}

// Create a formatted Djinn PG line to add to SAM header
func NewPG(hdr *sam.Header, cl string) *sam.Program {
	progs := hdr.Progs()

	referenced := make(map[string]bool, len(progs))
	seen := make(map[string]bool, len(progs))
	for _, p := range progs {
		if pp := p.Previous(); pp != "" {
			referenced[pp] = true
		}
		seen[p.UID()] = true
	}

	var prev string
	for _, p := range progs {
		if !referenced[p.UID()] {
			prev = p.UID()
			break
		}
	}

	uid := "djinn"
	for i := 1; seen[uid]; i++ {
		uid = fmt.Sprintf("djinn.%d", i)
	}

	return sam.NewProgram(uid, "djinn", cl, prev, "3.0")
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
	var w AlignmentWriter
	outChan := make(chan *sam.Record, cp)
	doneChan := make(chan bool)

	fh, err := os.Stdout, error(nil)
	if outFile != "-" {
		fh, err = os.Create(outFile)
		checkError(err)
	}

	bio := bufio.NewWriterSize(fh, buff)

	if uncompressed {
		sw, err := sam.NewWriter(bio, head, sam.FlagDecimal)
		checkError(err)
		w = &SamWriter{Writer: sw}
	} else {
		bw, err := bam.NewWriterLevel(bio, head, 4, threads)
		checkError(err)
		w = bw
	}

	go func() {
		for rec := range outChan {
			checkError(w.Write(rec))
		}
		w.Close()
		checkError(bio.Flush()) // always flush
		if fh != os.Stdout {
			fh.Close()
		}
		doneChan <- true
	}()

	return outChan, doneChan
}
