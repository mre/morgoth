package ymtn

import (
	"fmt"
	"github.com/golang/glog"
	"github.com/nvcook42/linalg"
	"github.com/nvcook42/linalg/lapack"
	"github.com/nvcook42/matrix"
	"github.com/nvcook42/morgoth/stat"
	"math"
)

type Event struct {
	Series string
	Data   []float64

	//Stats of the event data
	diffs []float64
	mean  float64
	std   float64

	//Flag indicating whether the stats have been computed for this event
	stats bool
}

const pearsonThreshold = 0.80
const epsilon = 1e-20

var nrhsOpt = linalg.IntOpt("nrhs", 1)
var ldBOpt = linalg.IntOpt("ldB", 1)

//Convert a set of motifs into events
// reducing the total set as necessary
func MotifsToEvents(motifs [][]Motif, T int, seriesName string) []Event {

	allEvents := make([]Event, 0, len(motifs))

	for _, pattern := range motifs {
		occurences := len(pattern)
		events := make([]Event, occurences)
		for i, motif := range pattern {
			events[i] = createEventFromMotifOccurence(motif, T, seriesName)
		}

		allEvents = append(allEvents, events...)
	}

	isChanged := true
	for isChanged {
		allEvents, isChanged = reduceOccurences(allEvents)
	}

	return allEvents
}

func createEventFromMotifOccurence(m Motif, T int, seriesName string) Event {
	data := make([]float64, T)
	for t := m.Beg; t <= m.End; t++ {
		data[t] = 1.0
	}
	event := Event{
		Series: seriesName,
		Data:   data,
	}
	return event
}

//Calculates the statistics needed to calculate the
// pearson correlation coefficient on event
func calcPearsonStats(event *Event) {
	count := float64(len(event.Data))
	sum := 0.0
	for _, x := range event.Data {
		sum += x
	}
	mean := sum / count
	variance := 0.0
	diffs := make([]float64, len(event.Data))
	for i, x := range event.Data {
		diff := x - mean
		diffs[i] = diff
		variance += diff * diff
	}

	event.diffs = diffs
	event.mean = mean
	event.std = math.Sqrt(variance)
	event.stats = true
}

//Calculate the Pearson correlation coefficient
// between two series.
// r = p(e1(t), e2(t+lag))
func pearsonCorrelation(e1, e2 *Event, lag int) float64 {

	if !e1.stats {
		calcPearsonStats(e1)
	}

	if !e2.stats {
		calcPearsonStats(e2)
	}

	l := len(e1.Data)

	sum := 0.0
	for i := 0; i < l; i++ {
		j := i + lag
		if j < 0 || j >= l {
			continue
		}
		sum += e1.diffs[i] * e2.diffs[j]
	}

	r := sum / (e1.std * e2.std)

	return r
}

func reduceOccurences(events []Event) ([]Event, bool) {

	//Compress all event occurences that are the same
	occurences := len(events)
	finalEvents := make([]Event, 0, occurences)
	isChanged := false
	combined := make([]bool, occurences)
	for i := 0; i < occurences; i++ {
		e1 := events[i]
		for j := i + 1; j < occurences; j++ {
			if combined[i] || combined[j] {
				continue
			}
			e2 := events[j]
			r := pearsonCorrelation(&e1, &e2, 0)
			if r > pearsonThreshold {
				isChanged = true
				combined[i] = true
				combined[j] = true
				e := Event{Series: e1.Series}
				e.Data = make([]float64, len(e1.Data))
				for k := range e1.Data {
					e.Data[k] = e1.Data[k] * e2.Data[k]
				}
				finalEvents = append(finalEvents, e)
			}
		}
	}

	for i := 0; i < occurences; i++ {
		if !combined[i] {
			finalEvents = append(finalEvents, events[i])
		}
	}

	return finalEvents, isChanged
}

// Find the dead time between e1 and e2
func findDeadTime(e1, e2 *Event, maxDeadTime int) int {
	max := 0.0
	maxLag := 0
	for lag := 0; lag < maxDeadTime; lag++ {
		r := pearsonCorrelation(e1, e2, lag)
		if r > max {
			max = r
			maxLag = lag
		}
	}
	return maxLag
}

func findPairwiseDeadTimes(events []Event) [][]int {

	count := len(events)
	deadTimes := make([][]int, count)

	for i, e1 := range events {
		deadTimes[i] = make([]int, count)
		for j, e2 := range events {
			if i == j {
				continue
			}
			d := findDeadTime(&e1, &e2, 30)
			deadTimes[i][j] = d
		}
	}

	return deadTimes
}

type Node struct {
	Series   string
	Begin    int
	Edges    []Edge
	IsHidden bool
}

type Edge struct {
	Cause  *Node
	Effect *Node
	Weight int
}

func (self *Node) Name() string {
	if !self.IsHidden{
		return fmt.Sprintf("%s_%d", self.Series, self.Begin)
	} else {
		return self.Series
	}
}

//Perform the Granger Causality tests
// on a set of event signals. The result
// is a list of nodes that belong to a weighted
// directed graph
func GrangerCausality(events []Event) []*Node {

	m := len(events)
	if m == 0 {
		return []*Node{}
	}
	T := len(events[0].Data)

	deadTimes := findPairwiseDeadTimes(events)

	nodes := make([]*Node, m)
	hiddenCount := 0

	for i := range nodes {
		nodes[i] = new(Node)
		node := nodes[i]

		node.Series = events[i].Series
		node.Edges = make([]Edge, 0)
		j := 0
		for ; events[i].Data[j] == 0 && j < T; j++ {
		}
		node.Begin = j
	}

	glog.V(2).Infoln("Nodes: ", nodes)

	for i := range events {
		for j := range events {
			if i == j {
				continue
			}
			determineCausality(i, j, events, deadTimes, nodes, &hiddenCount)
		}
	}

	return nodes
}

func determineCausality(i, j int, events []Event, deadTimes [][]int, nodes []*Node, hiddenCount *int) {

	m := len(events)
	p := 5

	T := len(events[j].Data)

	x0 := matrix.FloatZeros(T, 1)
	xlag0 := matrix.FloatZeros(T, (m-1)*p+1)
	x1 := matrix.FloatZeros(T, 1)
	xlag1 := matrix.FloatZeros(T, m*p+1)
	for t := 0; t < T; t++ {
		xlag1.SetAt(t, 0, 1.0)
		x1.SetAt(t, 0, events[j].Data[t])
		for l := 0; l < m; l++ {
			for k := 0; k < p; k++ {
				dt := deadTimes[l][j]
				tl := t - k - dt
				if tl >= 0 && tl < T {
					xlag1.SetAt(t, l*p+k+1, events[l].Data[tl])
				}
			}
		}

		xlag0.SetAt(t, 0, 1.0)
		x0.SetAt(t, 0, events[j].Data[t])
		for l := 0; l < m; l++ {
			if l == i {
				continue
			}
			for k := 0; k < p; k++ {
				dt := deadTimes[l][j]
				tl := t - k - dt
				if tl >= 0 && tl < T {
					adjL := l
					if l > i {
						adjL--
					}
					xlag0.SetAt(t, adjL*p+k+1, events[l].Data[tl])
				}
			}
		}

	}

	ssr1 := solveAutoRegression(x1, xlag1)
	glog.V(4).Infoln("SSR1", ssr1)

	ssr0 := solveAutoRegression(x0, xlag0)
	glog.V(4).Infoln("SSR0", ssr0)

	Sp, crit := calcTestStatAndCrit(ssr0, ssr1, float64(T), float64(p))

	glog.V(4).Infoln("Sp", Sp)
	glog.V(4).Infoln("crit", crit)

	if Sp > crit {
		x2 := matrix.FloatZeros(T, 1)
		xlag2 := matrix.FloatZeros(T, (m-1)*p+1+p)
		for t := 0; t < T; t++ {
			xlag2.SetAt(t, 0, 1.0)
			x2.SetAt(t, 0, events[j].Data[t])
			for l := 0; l < m; l++ {
				if l == i {
					continue
				}
				for k := 0; k < p; k++ {
					dt := deadTimes[l][j]
					tl := t - k - dt
					if tl >= 0 && tl < T {
						adjL := l
						if l > i {
							adjL--
						}
						xlag2.SetAt(t, adjL*p+k+1, events[l].Data[tl])
					}
				}
			}

			for k := 0; k < p; k++ {
				dt := deadTimes[j][i]
				tl := t + k + dt
				if tl >= 0 && tl < T {
					xlag2.SetAt(t, (m-1)*p+1+k, events[i].Data[tl])
				}
			}
		}

		ssr2 := solveAutoRegression(x2, xlag2)
		glog.V(4).Infoln("SSR2", ssr2)

		Sfp, crit2 := calcTestStatAndCrit(ssr0, ssr2, float64(T), float64(p))
		glog.V(4).Infoln("Sfp", Sfp)
		glog.V(4).Infoln("crit2", crit2)

		if Sfp < crit2 {
			// i -> j
			glog.V(2).Infoln("i -> j")
			addEdge(nodes[i], nodes[j], deadTimes[i][j])
		} else {
			// h -> i and h -> j
			glog.V(2).Infoln("h -> i and h -> j")
			hidden := Node{
				Series: fmt.Sprintf("hidden%d", *hiddenCount),
				IsHidden: true,
			}
			addEdge(&hidden, nodes[i], 0)
			addEdge(&hidden, nodes[j], 0)
			*hiddenCount = *hiddenCount + 1
		}
	}
	//Do nothing the relationship is not causal
}

func addEdge(cause, effect *Node, weight int) {
	edge := Edge{
		Cause:  cause,
		Effect: effect,
		Weight: weight,
	}
	cause.Edges = append(cause.Edges, edge)
}

func calcTestStatAndCrit(a, b, T, p float64) (float64, float64) {
	if b == 0 {
		return 0, 0
	}
	d1 := p
	d2 := T - 2*p - 1
	statistic := ((a - b) / d1) / (b / d2)
	crit := stat.F_CDF_At(d1, d2, 0.05)
	return statistic, crit
}

func solveAutoRegression(x, xlag *matrix.FloatMatrix) float64 {

	sn := int(math.Min(float64(xlag.Rows()), float64(xlag.Cols())))
	singularValues := make([]float64, sn)
	rank, err := lapack.Gelsd(xlag, x, singularValues, nrhsOpt)
	if err != nil {
		glog.Errorln(err)
	}
	glog.V(4).Infoln("Rank:", rank)

	ssr := 0.0
	for i := 0; i < x.Rows(); i++ {
		r := x.GetAt(i, 0)
		r2 := r * r
		if r2 > epsilon {
			ssr += r2
		}
	}

	return ssr
}
