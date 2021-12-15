### General information

This is a simple event generator dedicated to dilepton production process in ultra-peripheral collisions.

### Requirements

* Pythia8 event generator (for decays). See [Pythia website](https://pythia.org/) for installation instructions. It is
also possible to use Pythia6. Note that in both cases one needs to build ROOT with Pythia support. Note that Pythia is
only used for tau decays and completely optional
* ROOT (for calculations and Lorentz vectors). See [ROOT website](https://root.cern.ch/).
* CMAKE 2.8 (or newer) & gcc 4.8.5 (or newer)
* HepMC: It is possible to use HepMC output format. Source code and installation instructions 
can be obtained from [HepMC repository](https://gitlab.cern.ch/hepmc/HepMC3).
* Optionally: a compiler supporting OpenMP 4.5+ (some pragmas may be different for different versions).

### Other 3rd party libraries

This generator utilizes simple [plog](https://github.com/SergiusTheBest/plog) library for logging.

### Quick guide

* Clone this repository
* Install ROOT with Pythia8 support (recommended for decays) or Pythia6 support
* Load ROOT environment
* Setup Pythia8 environment by exporting `PYTHIA8=/path/to/pythia/install/dir`
* Optionally: build the generator with Pythia6 using cmake flag `BUILD_WITH_PYTHIA6=ON`
* Optionally: build the generator with Pythia8 using cmake flag `BUILD_WITH_PYTHIA8=ON`
* Optionally: build the generator with OpenMP support for parallel computation of two-photon luminosity
  (the most CPU-consuming operation) using cmake flag `BUILD_WITH_OPENMP=ON`
* Optionally: install HepMC library following instructions from the GitLab repository,
setup environment with `export HEPMC_ROOT=/path/to/hepmc/install/dir` and build the generator 
with HepMC support using `BUILD_WITH_HEPMC=ON`
* For details, see `CMakeLists.txt`
* Build the project:
```shell
cd path/to/cloned/repo
mkdir build
cd build
cmake ..
make
```

To run the generator use

```shell
./generate debug_level [optionally]number_of_threads
```

here, `debug_level` is `0`, `1`, or `2`. In debug mode, the generator will also print number of the event being
processed and verbose information about produced particles. In the most verbose mode (`2`) the program will also print
intermediate calculation results, so it is recommended to use it very carefully (and also to look in the code, if
possible).

The default number of threads is 1.

Generated events will be stored to `events_{a_lep}_{pt_cut}.root`, where `a_lep` is a value of the anomalous magnetic
moment and `pt_cut` is a minimal transverse momentum for a pair of leptons. If the generator is built with HepMC flag,
events will be written to `events_{a_lep}_{pt_cut}.hepmc2`

### Input parameters

Input parameters can be specified in the `parameters.in` file. The example can be found in the top directory of this
repository. Available parameters are the following:

```
NUCLEUS_Z 82      # atomic number of the incoming nuclei
NUCLEUS_A 208     # atomic mass of the incoming nuclei
WS_R 6.68         # Woods-Saxon parameters: R
WS_A 0.447        #                         a
SQRTS 5020        # sqrt(s) in the CM frame
LEP_PDG 15        # lepton code according to Monte Carlo numbering scheme from PDG
LEP_A 0           # lepton anomalous magnetic moment
NEVENTS 1000      # number of events to be generated
DO_PT_CUT 0       # enable pt cut: 0 -- off, 1 -- on
PT_MIN 0          # pt cut
ZMIN -1           # min. z = cos(theta) for the elementary cross section
ZMAX 1            # max. z = cos(theta) for the elementary cross section
MMIN 3.56         # min. m for the elementary/nuclear cross section
MMAX 50           # max. m for the elementary/nuclear cross section
YMIN -6           # min. y for the nuclear cross section
YMAX 6            # max. y for the nuclear cross section
BINS_Z 100        # cross section binnings: bins for z
BINS_M 1001       #                         bins for m
BINS_Y 121        #                         bins for y
FLUX_POINT 1      # use point flux approximation or not
PYTHIA_VERSION 8  # Pythia version: 6, 8 or -1. -1 means that Pythia will not be used at all
PYTHIA8_FSR 0     # For Pythia8 only: simulate final state radiation (EM showers)
PYTHIA8_DECAYS 0  # For Pythia8 only: switch to turn on/off lepton decays
SEED -1           # Seed for random numbers generator. '-1' -> random seed
```

* Note that the order, and the number of the parameters are not fixed.
* If a parameter is not specified by user, the default value will be used.
* IMPORTANT: make sure that you have built the generator with support of a desired Pythia version,
that you are going to pass via `PYTHIA_VERSION`
* Note that `PYTHIA8_FSR` and `PYTHIA8_DECAYS` only work for Pythia8. Decays are always enabled if
the generator runs with Pythia6.

### Tips

* The generator calculates two-photon luminosity and caches it into `hD2LDMDY.root`. 
This file will be picked automatically if found. The calculation process may take a lot of time,
so you may want to keep pre-calculated grid for further usage. Note that you need to recalculate
it in case if you have changed grid input parameters (e.g., binning and/or range in M/Y).