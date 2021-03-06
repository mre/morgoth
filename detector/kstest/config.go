package kstest

import (
	"github.com/nvcook42/morgoth/Godeps/_workspace/src/github.com/golang/glog"
	"github.com/nvcook42/morgoth/Godeps/_workspace/src/gopkg.in/validator.v2"
	"github.com/nvcook42/morgoth/defaults"
)

type KSTestConf struct {
	Confidence      uint `yaml:"confidence" validate:"min=1,max=5" default:"1"`
	NormalCount     uint `yaml:"normal_count" validate:"nonzero" default:"3"`
	MaxFingerprints uint `yaml:"max_fingerprints" validate:"nonzero" default:"20"`
}

func (self *KSTestConf) Default() {
	err := self.Validate()
	if err != nil {
		errs := err.(validator.ErrorMap)
		for fieldName := range errs {
			glog.Infof("Using default for KSTestConf.%s", fieldName)
			defaults.SetDefault(self, fieldName)
		}
	}

}

func (self *KSTestConf) Validate() error {
	return validator.Validate(self)
}
