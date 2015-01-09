package main

import (
	"github.com/golang/glog"
	"flag"
	"bufio"
	"os"
	"strings"
	"strconv"
	"github.com/nvcook42/morgoth/learners/ymtn"
)

var inFile = flag.String("in", "", "Input file. The first row should contain names of the series. Each successive row should contain numerical values. The file can have as many columns (time series) as desired within reason. All values should be space separated.")
var outFile = flag.String("out", "mnnet.dot", "Output file in the graphiz .dot format. The format is a simple text based format so it can be inpsected even without the graphiz program")



func main() {
	flag.Parse()

	if len(*inFile) == 0 {
		flag.Usage()
		return
	}

	infile, err := os.Open(*inFile)
	if err != nil {
		glog.Fatal(err)
	}
	defer infile.Close()

	scanner := bufio.NewScanner(infile)

	series, size := parseInFile(scanner)
	glog.Infof("Found %d series\n", len(series))

	allEvents := make([]ymtn.Event, 0)
	for s := range series {
		glog.Infoln("Preprocessing ", s)
		if len(series[s]) != size {
			glog.Fatalf("Series %s has %d data points not %d the expected number", s, len(series[s]), size)
		}
		scores := ymtn.RSST(series[s], 5, 4)
		motifs := ymtn.DetectMotifs(series[s], scores, 5, 50)
		events := ymtn.MotifsToEvents(motifs, size, s)
		allEvents = append(allEvents, events...)

	}

	glog.Infof("Found %d events\n", len(allEvents))

	wdg := ymtn.GrangerCausality(allEvents)

	glog.Infoln("Generating output of graph")
	outfile, err := os.Create(*outFile)
	if err != nil {
		glog.Fatal(err)
	}
	defer outfile.Close()

	out := bufio.NewWriter(outfile)
	ymtn.OutputGraph(out, wdg)
}

func parseInFile(in *bufio.Scanner) (map[string][]float64, int) {
	series := make(map[string][]float64)
	headers := make([]string, 0)
	rows := 0
	for in.Scan() {
		rows++
		line := in.Text()
		//Parse header row
		if len(headers) == 0 {
			names := strings.Split(line, " ")
			for _, name := range names {
				headers = append(headers, name)
				series[name] = make([]float64, 0, 100)
			}
			continue
		}

		//Parse data row
		values := strings.Split(line, " ")
		for i, value := range values {
			v,err := strconv.ParseFloat(value, 64)
			if err != nil {
				glog.Fatal(err)
			}
			series[headers[i]] = append(series[headers[i]], v)
		}
	}

	if err := in.Err(); err != nil {
		glog.Fatal(err)
	}

	return series, rows -1
}
