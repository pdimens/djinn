package xam

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// AlignmentReader is satisfied by both bam.Reader and sam.Reader
type AlignmentReader interface {
	Read() (*sam.Record, error)
	Header() *sam.Header
	Close() error
}

// wrapper for SAM to satisfy the interface
type SamReader struct {
	*sam.Reader
	closer io.Closer
}

func (s *SamReader) Close() error {
	return s.closer.Close()
}

// An interface that satisfies both SAM and BAM writers
type AlignmentWriter interface {
	Write(*sam.Record) error
	Close() error
}

// A wrapper to implement Close() to satsify the SAM/BAM writer interface.
type SamWriter struct {
	*sam.Writer
}

func (s *SamWriter) Close() error {
	return nil // sam.Writer has no resources to release
}

// Given an infile, will check for the gzip magic bytes and if present, will open an bam.Reader,
// otherwise opening a sam.Reader. Returns a reader that satisfied the AlignmentReader interface,
// a closer function to use with defer after invocation, and an error.
func OpenXAM(infile string, threads, buffer int) (AlignmentReader, func(), error) {
	var (
		bf  io.ReadCloser
		err error
	)

	switch infile {
	case "", "-":
		bf = os.Stdin
	default:
		bf, err = os.Open(infile)
		if err != nil {
			return nil, nil, fmt.Errorf("opening file: %w", err)
		}
	}

	bufReader := bufio.NewReaderSize(bf, buffer)

	magic, err := bufReader.Peek(2)
	if err != nil {
		if bf != os.Stdin {
			bf.Close()
		}
		return nil, nil, fmt.Errorf("peeking magic bytes for SAM/BAM determination: %w", err)
	}

	baseCleanup := func() {
		if bf != os.Stdin {
			bf.Close()
		}
	}

	if magic[0] == 0x1f && magic[1] == 0x8b {
		br, err := bam.NewReader(bufReader, threads)
		if err != nil {
			baseCleanup()
			return nil, nil, fmt.Errorf("creating BAM reader: %w", err)
		}
		return br, func() { br.Close(); baseCleanup() }, nil
	}

	sr, err := sam.NewReader(bufReader)
	if err != nil {
		baseCleanup()
		return nil, nil, fmt.Errorf("creating SAM reader: %w", err)
	}
	wrapped := &SamReader{Reader: sr, closer: bf}
	return wrapped, func() { wrapped.Close() }, nil
}
