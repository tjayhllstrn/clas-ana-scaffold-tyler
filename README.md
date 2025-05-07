## Installation

1. Clone the repository locally.

```
git clone https://github.com/Gregtom3/clas-ana-scaffold
```

2. Create a new python virtual environment (recommended: home directory).

```
python3 -m venv ~/.clas-ana-scaffold
```

3. Source the new environment

```
source ~/.clas-ana-scaffold/bin/activate.csh
```

or (depending on your shell)

```
source ~/.clas-ana-scaffold/bin/activate
```

4. Install python packages for this project

```
pip install -r requirements.txt
```

---

## Repository Structure

The structure of the repository is currently minimal. Assuming you are on `ifarm`, symoblic links such as `pass1_fall2018_rga_inbending_nSidis/` point to the location of `.hipo` files from a run period. 

In `src/`, C++ code handles some of the heavy-lifting involved in parsing/analyzing the `.hipo` file data. For example, in `Kinematics.C`, several methods are written to calculate DIS, single hadron, and dihadron kinematics. In `CutManager.C`, an array of particle/event cuts are defined, similar to those used in the RG-A analysis note. In the subdirectory `src/gbt/`, python code responsible for loading and running the Gradient Boosted Trees algorithm on photons is included. The models presented are pre-trained by Gregory and should be reliable enough for `RG-A`, `RG-B`, and `RG-C` (so long as the beam energy is in the ballpark of 10.6 GeV).

In `macros/`, there are pairs of `hipo2tree` and `Builder` codes. For now, only a `pi0` (single hadron SIDIS) code is included. The `hipo2tree` code handles reading the `.hipo` files, applying particle/event cuts, and storing the results in a TTree called `EventTree`. The `EventTree`'s rows correspond to separate events. The TBranches store either single values, such as the Bjorken x of the event, or arrays, such as the x-momentum of each particle saved in that event. In Pi0-related analyses, a followup TTree called `MLinput` is created, which serves as the input to the machine learning model for predicting if the photons are true/false.

The `Builder` code then parses event-by-event data in `EventTree` to reconstruct single hadron/dihadron kinematics. In the case of Pi0-related studies, an intermediate step calling the Gradient Boosted Trees algorithm is needed. This step creates a new per-particle branch called `p_gamma` in the `EventTree` which the `Builder` code may use to place further cuts.

In `scripts/`, Ruby scripts are provided to run the pipeline all in one go. This way you can run a single line of code and generate TTree's containing analysis-worthy DIS, SIDIS, or DISIDIS kinematics. An example is provided below.

**NOTE: YOU MUST BE SOURCING A PROPER PYTHON ENVIRONMENT FOR THE GBT ALGORITHM TO HAVE THE PROPER INCLUDES. SEE INSTALLATION FOR A WALKTHROUGH.**

## Example A (Pi0 Analysis Pipeline)

```sh
./scripts/run_single_hipo.rb -t pi0 -i pass1_fall2018_rga_inbending_nSidis/nSidis_005032.hipo -o out/test -n 500
```

* `-t, --type=TYPE`       : SIDIS type (only `pi0` supported)
* `-i, --input=FILE`      : Path to the input `.hipo` file
* `-o, --outdir=DIR`      : Directory where output `.root` and results will be saved
* `-n, --max-events=N`    : Maximum number of events to process (default: `100`, use `-1` for all events)
* `-h, --help`            : Show help message

This will:

1. Create `out/test/` if it doesn't exist.
2. Creates event-by-event TTree from `nSidis_005032.hipo` → `out/test/nSidis_005032.root`.
3. Run GBT photon-predictions using the appropriate model.
4. Execute `pi0Builder` on `out/test/nSidis_005032.root`.

