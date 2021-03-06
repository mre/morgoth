package mocks

import "github.com/nvcook42/morgoth/engine"
import "github.com/nvcook42/morgoth/Godeps/_workspace/src/github.com/stretchr/testify/mock"

import "github.com/nvcook42/morgoth/schedule"

type Engine struct {
	mock.Mock
}

func (m *Engine) Initialize() error {
	ret := m.Called()

	r0 := ret.Error(0)

	return r0
}
func (m *Engine) ConfigureSchedule(schedule *schedule.Schedule) error {
	ret := m.Called(schedule)

	r0 := ret.Error(0)

	return r0
}
func (m *Engine) GetReader() engine.Reader {
	ret := m.Called()

	r0 := ret.Get(0).(engine.Reader)

	return r0
}
func (m *Engine) GetWriter() engine.Writer {
	ret := m.Called()

	r0 := ret.Get(0).(engine.Writer)

	return r0
}
