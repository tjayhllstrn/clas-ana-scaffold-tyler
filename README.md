1. **HIPO → ROOT TTree** using the `hipo2tree_pi0.C` macro (via `clas12root`).
2. **Gradient-boosted tree (GBT) photon-prediction** with `src/gbt/predict.py`. This step creates a new branch `p_gamma` for each photon.
3. **π⁰ building** using the `pi0Builder.C` macro on the resulting ROOT file.

---

## Usage

```sh
./run_single_hipo.rb -t pi0 -i path/to/input.hipo -o path/to/output/dir [options]
```

### Options

* `-t, --type=TYPE`       : SIDIS type (only `pi0` supported)
* `-i, --input=FILE`      : Path to the input `.hipo` file
* `-o, --outdir=DIR`      : Directory where output `.root` and results will be saved
* `-n, --max-events=N`    : Maximum number of events to process (default: `100`, use `-1` for all events)
* `-h, --help`            : Show help message

## Example

```sh
./run_single_hipo.rb -t pi0 -i data/nSidis_005032.hipo -o out/test -n 500
```

This will:

1. Create `out/test/` if it doesn't exist.
2. Creates event-by-event TTree from `nSidis_005032.hipo` → `out/test/nSidis_005032.root`.
3. Run GBT photon-predictions using the appropriate model.
4. Execute `pi0Builder` on `out/test/nSidis_005032.root`.

