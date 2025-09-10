package fastqreader

import (
	"iter"
	"log"
	"math/rand"
	"regexp"
)

var standardRe = regexp.MustCompile(`BX:Z:(\S+)\s`)
var vxRe = regexp.MustCompile(`VX:i:([01])\s`)
var invalidRe = regexp.MustCompile(`^0_|_0_|_0$|00|N`)
var haplotaggingRe = regexp.MustCompile(`BX:Z:(A\d{2}C\d{2}B\d{2}D\d{2})\s`)
var stlfrRe = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)\s`)
var tellseqRe = regexp.MustCompile(`:([ATCGN]+)\s`)

/*
Parse through the first 200 records of a paired-end fastq and figure out what kind of
format it's in. Returns early if a format is detected.
*/
func FindFastqFormat(r1, r2 string) string {
	var err error
	var barcode_type string = "unknown"

	fqr, err := OpenFastQ(r1, r2)

	if err != nil {
		log.Fatal(err)
	}

	for range 200 {
		barcode_type, err = fqr.InferReadType()
		if err != nil {
			log.Fatal("Error processing the input FASTQ files.")
		}
		if barcode_type != "unknown" {
			return barcode_type
		}
	}
	log.Fatal("Unable to identify the input FASTQ files as one of haplotagging, stLFR, or TELLseq. Unable to proceed.")
	return "unknown"
}

// Product generates the Cartesian product of input slices
func product[T any](inputs ...[]T) iter.Seq[[]T] {
	return func(yield func([]T) bool) {
		if len(inputs) == 0 {
			return
		}

		// Check if any input is empty
		for _, input := range inputs {
			if len(input) == 0 {
				return
			}
		}

		indices := make([]int, len(inputs))
		result := make([]T, len(inputs))

		for {
			// Build current combination
			for i, idx := range indices {
				result[i] = inputs[i][idx]
			}

			// Yield current combination
			if !yield(result) {
				return
			}

			// Generate next indices (like odometer)
			carry := 1
			for i := len(indices) - 1; i >= 0 && carry > 0; i-- {
				indices[i] += carry
				if indices[i] >= len(inputs[i]) {
					indices[i] = 0
					carry = 1
				} else {
					carry = 0
				}
			}

			// If we've cycled through all combinations
			if carry > 0 {
				break
			}
		}
	}
}

// Create a slice of N repetitions of value V
func sliceRepeat[T []string](N int, V T) []T {
	retval := make([]T, 0, N)
	for range N {
		Vshuff := V
		for i := len(Vshuff) - 1; i > 0; i-- {
			j := rand.Intn(i + 1)
			Vshuff[i], Vshuff[j] = Vshuff[j], Vshuff[i]
		}
		retval = append(retval, Vshuff)
	}
	return retval
}

// StringProduct is a convenience function for string combinations
func BarcodeGenerator() iter.Seq[string] {
	return func(yield func(string) bool) {
		bases := []string{"A", "T", "C", "G"}
		inputs := sliceRepeat(20, bases)
		for combo := range product(inputs...) {
			result := ""
			for _, s := range combo {
				result += s
			}
			if !yield(result) {
				return
			}
		}
	}
}
