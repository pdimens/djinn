package count

import (
	"bufio"
	"fmt"
	"os"
	"strconv"

	"djinn/xam"
)

func CountXam(infile string, invalid bool, threads int) error {
	infile = xam.FileOrStdin(infile)
	set := make(map[string]int16, 7_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)

	// ── open writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if !hasBX {
			continue
		}
		if !vxVal && !invalid {
			continue
		}
		set[bxVal]++
	}

	for key, val := range set {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.WriteString(strconv.Itoa(int(val)))
		writer.WriteByte('\n')
		fmt.Printf("%s\t%d\n", key, val)
	}
	return nil
}
