# How to run

Everything below is tested on Ubuntu 24.04 with Python 3.10–3.13 and by the CI workflow
([`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

## Requirements

- a C++17 compiler (GCC ≥ 9 or Clang ≥ 10), CMake ≥ 3.20
- Eigen 3.4 (`sudo apt install libeigen3-dev`; if it is missing, CMake downloads the headers)
- Python ≥ 3.10 with `pip`
- network access on the first build: Catch2 v3 is fetched by CMake for the tests

## C++ library, tests and examples

```bash
make build        # cmake -S . -B build/cmake && cmake --build build/cmake
make test         # ctest: every algorithm against a closed form, a hand-worked example or brute force
```

The examples are small command-line programs, one per chapter:

```bash
./build/cmake/cpp/examples/01_nmea_summary data/gnss/nmea_log2.nmea 10
```

## Python module

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[notebooks]"             # builds the C++ core and the pybind11 module
pip install -e ".[notebooks,reference]"   # + Open3D, used as a reference in some notebooks
```

`pip install -e .` compiles the C++ code once; after changing C++ files run it again (`make install`).

```python
import sensor_fusion as sf
log = sf.nmea.parse_log_file(str(sf.data.path("gnss/nmea_log2.nmea")))
```

## Notebooks

Each chapter has two copies of its notebook:

- `2_notebooks/exercises/NN_*.ipynb` — the reader's copy, with ✏️ cells to fill in;
- `2_notebooks/solutions/NN_*.ipynb` — solved, executes top to bottom.

```bash
jupyter lab 2_notebooks/exercises      # work through a chapter
make notebooks                         # execute all solution notebooks (stores outputs)
make exercises                         # regenerate exercises/ from solutions/ (tools/make_exercises.py)
```

### Colab

Open a notebook from GitHub in Colab (`https://colab.research.google.com/github/Rudip1/learn-sensor-fusion/blob/main/2_notebooks/exercises/01_gnss_fundamentals.ipynb`).
The first cell installs the module with
`%pip install -q git+https://github.com/Rudip1/learn-sensor-fusion` (a few minutes: it compiles the C++).
Data files are downloaded on first use into `~/.cache/sensor_fusion` (`sensor_fusion.data.path`).

## Figures

The figures in `1_theory/figures/` are generated, never edited by hand:

```bash
make figures      # python tools/make_figures.py; prints the numbers quoted in the theory text
```

## Everything

```bash
make all          # build + test + install + notebooks
```
