package ymtn_test

import (
	"flag"
	"github.com/golang/glog"
	"github.com/nvcook42/linalg/lapack"
	"github.com/nvcook42/matrix"
	"github.com/nvcook42/morgoth/learners/ymtn"
	"github.com/stretchr/testify/assert"
	"testing"
)

func init() {
	flag.Parse()
	if testing.Verbose() {
		flag.Set("logtostderr", "true")
		flag.Set("vmodule", "events=2,events_test=3")
	}
}

func TestMotifToEvents(t *testing.T) {
	assert := assert.New(t)

	size := 100
	x := make([]float64, size)
	for i := range x {
		//Saw tooth pattern that has clear change points
		x[i] = float64((i + 1) % 10)
	}
	scores := ymtn.RSST(x, 5, 4)

	motifs := ymtn.DetectMotifs(x, scores, 5, 25)

	events := ymtn.MotifsToEvents(motifs, size, "x")
	assert.NotNil(events)
	glog.V(1).Infoln("#Events:", len(events))
	glog.V(1).Infoln("Events:", events)
}

func TestGranger(t *testing.T) {
	assert := assert.New(t)

	size := 200
	x := make([]float64, size)
	for i := range x {
		//Saw tooth pattern that has clear change points
		x[i] = float64((i + 1) % 10)
	}
	y := make([]float64, size)
	for i := range y {
		//Saw tooth pattern that has clear change points
		y[i] = float64((i + 1) % 40)
	}
	xScores := ymtn.RSST(x, 5, 4)
	xMotifs := ymtn.DetectMotifs(x, xScores, 5, 25)

	xEvents := ymtn.MotifsToEvents(xMotifs, size, "x")
	glog.V(1).Infoln("# X Events:", len(xEvents))

	yScores := ymtn.RSST(y, 5, 4)
	yMotifs := ymtn.DetectMotifs(y, yScores, 5, 25)

	yEvents := ymtn.MotifsToEvents(yMotifs, size, "y")
	glog.V(1).Infoln("# Y Events:", len(yEvents))

	events := append(xEvents, yEvents...)

	dag := ymtn.GrangerCausality(events)
	assert.NotNil(dag)
	glog.V(1).Infoln("DAG:", dag)

}

func TestGels(t *testing.T) {
	assert := assert.New(t)

	a := matrix.FloatNew(4, 5, []float64{
		0.12, -6.91, -3.33, 3.97,
		-8.19, 2.22, -8.94, 3.33,
		7.69, -5.12, -6.72, -2.74,
		-2.26, -9.08, -4.40, -7.92,
		-4.71, 9.96, -9.98, -3.20,
	})

	b := matrix.FloatNew(5, 3, []float64{
		7.30, 1.33, 2.68, -9.62, 0.00,
		0.47, 6.58, -1.71, -0.79, 0.00,
		-6.28, -3.42, 3.46, 0.41, 0.00,
	})

	glog.V(1).Infoln("A:", a)
	glog.V(1).Infoln("B:", b)
	s := make([]float64, 4)
	rank, err := lapack.Gelsd(a, b, s)
	if err != nil {
		glog.Errorln(err)
	}

	assert.Equal(4, rank)

	glog.V(1).Infoln("SolA:", a)
	glog.V(1).Infoln("Solx:", b)
	glog.V(1).Infoln("Sols:", s)

}
