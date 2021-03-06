package config_test

import (
	"github.com/nvcook42/morgoth/Godeps/_workspace/src/github.com/stretchr/testify/assert"
	"github.com/nvcook42/morgoth/config"
	"testing"
)

func TestConfigShouldNotParseInvalidYAML(t *testing.T) {
	assert := assert.New(t)
	var data = `
---
data_engine:
  mongodb:
	use_sharding: false
`
	config, err := config.Load([]byte(data))
	assert.NotNil(err)
	assert.Nil(config)
}

func TestConfigShouldNotValidateBadEngineConf(t *testing.T) {
	assert := assert.New(t)

	var data = `
---
data_engine:
  bad_key: 1
`
	_, err := config.Load([]byte(data))
	assert.NotNil(err)

}

func TestConfigShouldNotValidateBadMetricConf(t *testing.T) {
	assert := assert.New(t)

	var data = `
---
metrics:
  - {}
  - {}
`
	_, err := config.Load([]byte(data))
	assert.NotNil(err)

}

func TestConfigShouldNotValidateBadFittingConf(t *testing.T) {
	assert := assert.New(t)

	var data = `
---
fittings:
	unknown: {}
	bad: {}
`
	_, err := config.Load([]byte(data))
	assert.NotNil(err)

}
