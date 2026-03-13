package fqtools

import (
	"bufio"
	"io"
	"os"
	"strings"

	"github.com/klauspost/pgzip"
)

// ── FASTQ reader ──────────────────────────────────────────────────────────────

type fastqRecord struct {
	name    string
	comment string
	seq     []byte
	qual    []byte // ASCII Phred+33
}

type readPair struct {
	fq1, fq2 fastqRecord
}

type fastqReader struct{ r *bufio.Reader }

func newFastqReader(r io.Reader) *fastqReader {
	return &fastqReader{r: bufio.NewReaderSize(r, 1<<19)}
}

func (f *fastqReader) next() (fastqRecord, bool) {
	header, err := f.r.ReadString('\n')
	if err != nil {
		return fastqRecord{}, false
	}
	header = strings.TrimRight(header, "\r\n")
	if len(header) == 0 || header[0] != '@' {
		return fastqRecord{}, false
	}
	header = header[1:]

	var name, comment string
	if before, after, ok := strings.Cut(header, " "); ok {
		name, comment = before, after
	} else {
		name = header
	}

	seqLine, _ := f.r.ReadString('\n')
	seqLine = strings.TrimRight(seqLine, "\r\n")
	_, _ = f.r.ReadString('\n') // skip '+'
	qualLine, _ := f.r.ReadString('\n')
	qualLine = strings.TrimRight(qualLine, "\r\n")

	return fastqRecord{name: name, comment: comment, seq: []byte(seqLine), qual: []byte(qualLine)}, true
}

// ── file helpers ──────────────────────────────────────────────────────────────

func openFastq(path string, gzBlocks int) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(path, ".gz") {
		gz, err := pgzip.NewReaderN(f, 3<<20, gzBlocks)
		if err != nil {
			f.Close()
			return nil, err
		}
		return struct {
			io.Reader
			io.Closer
		}{gz, multiCloser{gz, f}}, nil
	}
	return f, nil
}

type multiCloser []io.Closer

func (mc multiCloser) Close() error {
	var last error
	for _, c := range mc {
		if err := c.Close(); err != nil {
			last = err
		}
	}
	return last
}

// ─── FASTQ writer ─────────────────────────────────────────────────────────────

// fastqWriter wraps a pgzip.Writer with a bufio layer.
// Writes go directly into the bufio buffer, which flushes to pgzip only
// when its 1MB capacity is reached — no intermediate per-record copy.
type fastqWriter struct {
	gz  *pgzip.Writer
	buf *bufio.Writer
	dir string
}

func newFastqWriter(path string, readnum string, level, threads int) (*fastqWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	gz, err := pgzip.NewWriterLevel(f, level)
	if err != nil {
		f.Close()
		return nil, err
	}
	gz.SetConcurrency(2<<20, threads)
	return &fastqWriter{
		gz:  gz,
		buf: bufio.NewWriterSize(gz, 2<<20),
		dir: readnum,
	}, nil
}

// writeRecord writes a single FASTQ record directly into the bufio buffer.
// qual is raw Phred scores (0-40) as stored in sam.Record — +33 applied inline.
// bufio only flushes to pgzip when its 1MB buffer is full, so no intermediate
// per-record copy is needed.
func (fw *fastqWriter) writeRecord(name string, bxTag string, vxTag int, seq []byte, qual []uint8) error {
	// Header: @name\tVX:i:<vx>\tBX:Z:<bx>
	fw.buf.WriteByte('@')
	fw.buf.WriteString(name)
	fw.buf.WriteString(fw.dir)
	fw.buf.WriteString("\tVX:i:")
	writeInt(fw.buf, vxTag)
	fw.buf.WriteString("\tBX:Z:")
	fw.buf.WriteString(bxTag)
	fw.buf.WriteByte('\n')

	// Sequence
	fw.buf.Write(seq)
	fw.buf.WriteString("\n+\n")

	// Quality — convert raw Phred → ASCII Phred+33 inline
	for _, q := range qual {
		fw.buf.WriteByte(q + 33)
	}
	fw.buf.WriteByte('\n')
	return nil
}

func (fw *fastqWriter) close() error {
	if err := fw.buf.Flush(); err != nil {
		return err
	}
	return fw.gz.Close()
}

// writeInt writes a non-negative integer to b without allocating.
func writeInt(b *bufio.Writer, n int) {
	if n == 0 {
		b.WriteByte('0')
		return
	}
	var tmp [10]byte
	i := len(tmp)
	for n > 0 {
		i--
		tmp[i] = byte('0' + n%10)
		n /= 10
	}
	b.Write(tmp[i:])
}
