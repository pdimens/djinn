// Standalone benchmark driver for sort.SortByBX. Not part of the djinn
// build; run with `go run`. Reports wall time and peak RSS (via
// /proc/self/status's VmHWM, read just before exit) for a single sort run
// in its own process.
package main

import (
	"bufio"
	"djinn/sort"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"time"
)

func peakRSSKB() int64 {
	f, err := os.Open("/proc/self/status")
	if err != nil {
		return -1
	}
	defer f.Close()
	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		if strings.HasPrefix(line, "VmHWM:") {
			fields := strings.Fields(line)
			if len(fields) >= 2 {
				v, _ := strconv.ParseInt(fields[1], 10, 64)
				return v
			}
		}
	}
	return -1
}

func main() {
	in := flag.String("in", "", "input BAM path")
	out := flag.String("out", "", "output BAM path")
	tmp := flag.String("tmp", "", "temp dir")
	threads := flag.Int("threads", 4, "worker threads")
	label := flag.String("label", "", "label for this run, printed in the report line")
	flag.Parse()

	if *in == "" || *out == "" || *tmp == "" {
		log.Fatal("usage: sortrun -in=... -out=... -tmp=... [-threads=N] [-label=...]")
	}

	start := time.Now()
	if err := sort.SortByBX(*in, *out, *tmp, *threads, false); err != nil {
		log.Fatalf("SortByBX failed: %v", err)
	}
	elapsed := time.Since(start)

	rssKB := peakRSSKB()
	fmt.Printf("label=%s threads=%d elapsed_s=%.3f peak_rss_mb=%.1f\n",
		*label, *threads, elapsed.Seconds(), float64(rssKB)/1024)
}
