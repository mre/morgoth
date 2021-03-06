package generator

import (
	"github.com/nvcook42/morgoth/Godeps/_workspace/src/github.com/golang/glog"
	"github.com/nvcook42/morgoth/engine"
	metric "github.com/nvcook42/morgoth/metric/types"
	"github.com/nvcook42/morgoth/schedule"
	"math"
	"time"
)

var (
	TZero        = time.Unix(0, 0)
	tZeroSeconds = TZero.Unix()
)

type Ft func(t int64) float64

type GeneratorEngine struct {
	functions map[metric.MetricID]Ft
}

func New(functions map[metric.MetricID]Ft) *GeneratorEngine {
	ge := GeneratorEngine{
		functions,
	}

	return &ge
}

func (self *GeneratorEngine) Initialize() error {
	return nil
}

func (self *GeneratorEngine) GetReader() engine.Reader {
	return self
}

func (self *GeneratorEngine) GetWriter() engine.Writer {
	return self
}

func (self *GeneratorEngine) ConfigureSchedule(schedule *schedule.Schedule) error {
	return nil
}

//////////////////////
// Writer Methods
//////////////////////

func (self *GeneratorEngine) Insert(datetime time.Time, metric metric.MetricID, value float64) {
}

func (self *GeneratorEngine) RecordAnomalous(metric metric.MetricID, start, stop time.Time) {
}

func (self *GeneratorEngine) DeleteMetric(metric metric.MetricID) {
}

//////////////////////
// Reader Methods
//////////////////////

func (self *GeneratorEngine) GetMetrics() []metric.MetricID {
	return nil
}

func (self *GeneratorEngine) GetData(rotation *schedule.Rotation, metric metric.MetricID, start, stop time.Time) []engine.Point {
	data := make([]engine.Point, 0)

	f, ok := self.functions[metric]
	if !ok {
		glog.Warning("Queried data for undefined metric ", metric)
		return data
	}

	t := start
	for t.Before(stop) || t.Equal(stop) {
		value := f(t.Unix() - tZeroSeconds)
		data = append(data, engine.Point{t, value})
		t = t.Add(rotation.Resolution)
	}
	return data
}

func (self *GeneratorEngine) GetAnomalies(metric metric.MetricID, start, stop time.Time) []engine.Anomaly {
	return nil
}
func (self *GeneratorEngine) GetHistogram(rotation *schedule.Rotation, metric metric.MetricID, nbins uint, start, stop time.Time, min, max float64) *engine.Histogram {
	data := self.GetData(rotation, metric, start, stop)
	hist := new(engine.Histogram)
	hist.Count = uint(len(data))
	hist.Min = min
	hist.Max = max
	hist.Bins = make([]float64, nbins)

	step := (max - min) / float64(nbins)
	count := float64(hist.Count)

	for _, point := range data {
		v := point.Value
		if v < min || v > max {
			continue
		}
		i := int(math.Floor((v - min) / step))
		hist.Bins[i] += 1.0 / count
	}

	return hist
}
func (self *GeneratorEngine) GetPercentile(rotation *schedule.Rotation, metric metric.MetricID, percentile float64, start, stop time.Time) float64 {
	return 0.0
}
