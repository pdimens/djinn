package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
)

// check if the file exists and the directory is writable, return an
// error if not, otherwise nil
func filecheck(filename string) error {
	//var err error
	if err := fileExists(filename); err != nil {
		return err
	}
	if err := fileReadable(filename); err != nil {
		return err
	}
	return nil
}

// check if file exists. Returns a true if it does, otherwise false.
func fileExists(filename string) error {
	info, err := os.Stat(filename)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("\033[33;1m%s\033[0m does not exist", filename)
		}
		return err
	}
	if info.IsDir() {
		return fmt.Errorf("\033[33;1m%s\033[0m is a directory", filename)
	}
	return nil
}

// check if the input file has read permissions
func fileReadable(filename string) error {
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("\033[33;1m%s\033[0m does not have read persmissions.\n", filename)
	}
	file.Close()
	return nil
}

// check if executable 'e' is in path, returns true if it is, otherwise false
func checkIfExecInPath(e string) bool {
	// ignore output, only use error
	_, err := exec.LookPath(e)
	// return true if error nil, else false
	return err == nil
}

// Check if the input path has write permissions
func ensureWritableDir(path string) error {
	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.CreateTemp(dir, ".writetest")
	if err != nil {
		return err
	}
	f.Close()
	return os.Remove(f.Name())
}

// Parse the thread count provided to `threads` at the CLI and clamp it between 1 and the system max
func safethreads(threads int) int {
	maxCores := runtime.NumCPU()
	// clamp between 1 and max system threads
	return min(maxCores, max(threads, 1))
}
