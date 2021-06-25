### General information

This is a simple event generator dedicated to dilepton production process in ultra-peripheral collisions.

### Requirements

* Pythia8 event generator (for decays). See [Pythia website](https://pythia.org/) for installation instructions
* ROOT (for calculations and Lorentz vectors). See [ROOT website](https://root.cern.ch/).
* CMAKE 2.8 (or newer) & gcc 4.8.5 (or newer)

### Other 3rd party libraries

This generator utilizes simple [plog](https://github.com/SergiusTheBest/plog) library for logging.

### Quick guide

* Clone this repository, then build the project:

```shell
cd path/to/cloned/repo
mkdir build
cd build
cmake ..
make
```

To run the generator use

```shell
./generate debug_level
```

here, `debug_level` is `0` or `1` (for now). In debug mode, the generator will also print number of the event being processed and verbose information about produced particles.

Generated events will be stored to `events_{a_lep}_{pt_cut}.root`, where `a_lep` is a value of the anomalous magnetic moment and `pt_cut` is a minimal transverse momentum for a pair of leptons.

### Input parameters

Input parameters can be specified in the `parameters.in` file. The example can be found in the top directory of this
repository. Available parameters are the following:

```
NEVENTS 1000      # number of events to be generated
SQRTS 5020        # sqrt(s) in the CM frame
LEP_MASS 1.77682  # lepton mass
LEP_A 0           # lepton anomalous magnetic moment
DO_PT_CUT 0       # enable pt cut: 0 -- off, 1 -- on
PT_MIN 0          # pt cut
ZMIN -1           # min. z = cos(theta) for the elementary cross section
ZMAX 1            # max. z = cos(theta) for the elementary cross section
WMIN 3.56         # min. w for the elementary/nuclear cross section
WMAX 50           # max. w for the elementary/nuclear cross section
YMIN -6           # min. y for the nuclear cross section
YMAX 6            # max. y for the nuclear cross section
BINS_Z 100        # cross section binnings: bins for z
BINS_W 1001       #                         bins for w
BINS_Y 121        #                         bins for y
WS_RHO0 0.159538  # Woods-Saxon parameters: rho0
WS_RAD 6.68       #                         R
WS_A 0.447        #                         a
NUCLEUS_Z 82      # atomic number of the incoming nuclei
NUCLUES_A 208     # atomic mass of the incoming nuclei
```

* Note that the order, and the number of the parameters are not fixed.
* If a parameter is not specified by user, the default value will be used.

### Tips

* The generator calculates two-photon luminosity and caches it into `hD2LDMDY.root`. This file will be picked automatically if found. The calculation process may take a lot of time, so you may want to keep pre-calculated grid for further usage. Note that you need to recalculate it in case if you have changed grid input parameters (e.g., binning and/or range in z/w/Y).