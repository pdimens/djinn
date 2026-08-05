package extract

import (
	"bufio"
	"os"

	"djinn/xam"
)

func Extract(infile string, invalid bool, threads int) {
	infile = xam.FileOrStdin(infile)
	set := make(map[string]struct{}, 7_000_000)

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
		if _, ok := set[bxVal]; ok {
			continue
		} else {
			set[bxVal] = struct{}{}
			writer.WriteString(bxVal)
			writer.WriteByte('\n')
		}
	}
}
