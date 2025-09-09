package fastqreader

import (
	"io"
	"log"
	"strings"
)

/*Convert the forward-read part of a fastq record to a string and write it*/
func writeR1FastqRecord(record FastQRecord, barcode string, gzip_proc io.WriteCloser) {
	var fq_fmt string

	fq_fmt = record.Header + "/1\n"
	fq_fmt += barcode + strings.Repeat("N", 10) + string(record.Read1)
	fq_fmt += "\n+\n"
	fq_fmt += string(record.ReadQual1) + "\n"

	if _, err := gzip_proc.Write([]byte(fq_fmt)); err != nil {
		log.Fatalf("Error writing to gzip: %v", err)
	}
}

/*Convert the forward-read part of a fastq record to a string and write it*/
func writeR2FastqRecord(record FastQRecord, barcode string, gzip_proc io.WriteCloser) {
	var fq_fmt string

	fq_fmt = record.Header + "/2\n"
	fq_fmt += barcode + strings.Repeat("N", 10) + string(record.Read2)
	fq_fmt += "\n+\n"
	fq_fmt += string(record.ReadQual2) + "\n"

	if _, err := gzip_proc.Write([]byte(fq_fmt)); err != nil {
		log.Fatalf("Error writing to gzip: %v", err)
	}
}
