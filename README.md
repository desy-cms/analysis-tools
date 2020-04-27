# Analysis/Tools

**Core codes for the DESY CMS MSSM Higgs -> bb analysis framework**

* [Installation](#installation)
* [Calibrations](#calibrations)
* [Ntuples](#ntuples)
* [Example](#example)

## Installation

The codes as well as the ntuples are independent of CMSSW. However, in order to compile it uses `scram`.
So the latest version in any architecture should be fine. 

```bash
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src
cmsenv
git clone https://github.com/desy-cms/analysis-tools.git Analysis/Tools 
scram b -j4 USER_CXXFLAGS="-Wno-misleading-indentation"
```

:zzz: The `USER_CXXFLAGS="-Wno-misleading-indentation"` prevents a large number of warnings
from misleading identation in modules of the boost library. User can also define the
environment variable below in `.bashrc` or every time after the command `cmsenv`
```bash
export USER_CXXFLAGS="-Wno-misleading-indentation"
```

## Calibrations

Scale factors, efficiencies etc can be obtained from the [analysis-calibrations](https://github.com/desy-cms/analysis-calibrations) repository (see also the README.md in each run period directory). It is recommended to install it in the Analysis/Tools/data directory, for the framework loads the calibrations files from that path

```bash
cd $CMSSW_BASE/src/Analysis/Tools/data
git clone https://github.com/desy-cms/analysis-calibrations.git calibrations
```

## Ntuples

The lists of ntuples files can be obtained from the [analysis-ntuples](https://github.com/desy-cms/analysis-ntuples.git) repository (see also the README.md in each run period directory). The repository can be installed in a directory of your convenience, e.g.

```bash
cd $CMSSW_BASE/src/Analysis/Tools/test
git clone https://github.com/desy-cms/analysis-ntuples.git ntuples
```

## Example

A simple example macro can be found in [Analysis/Tools/bin/AnalyserSimpleExample.cc](bin/AnalyserSimpleExample.cc).
The example is a simple offline and online dijet selection (+ muon selection in the semileptonic case)
using signal MC samples and triggers from 2017 and 2018 periods. The macro uses a configuration file as
input. Configuration files are avaialble for both the all-hadronic and semileptonic cases for both 2017
and 2018 periods:
* [analyser_example_allhad_2017.cfg](test/analyser_example_allhad_2017.cfg)
* [analyser_example_allhad_2018.cfg](test/analyser_example_allhad_2018.cfg)
* [analyser_example_semilep_2017.cfg](test/analyser_example_semilep_2017.cfg)
* [analyser_example_semilep_2018.cfg](test/analyser_example_semilep_2018.cfg)

