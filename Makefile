# Build, test and run everything.   make all  =  build + test + install + notebooks
PYTHON    ?= python3
BUILD_DIR ?= build/cmake
JOBS      ?= $(shell nproc 2>/dev/null || echo 4)
PIP_FLAGS ?=

.PHONY: all build test install notebooks exercises check-exercises figures format clean

all: build test install notebooks

build:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD_DIR) -j $(JOBS)

test: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure -j $(JOBS)

# Python module (C++ core + bindings) in editable mode; re-run after changing C++ code.
install:
	$(PYTHON) -m pip install $(PIP_FLAGS) -e ".[notebooks]"

# Execute every solution notebook top to bottom and keep the outputs.
notebooks:
	$(PYTHON) tools/run_notebooks.py --inplace

# Regenerate the reader's notebooks from the solutions.
exercises:
	$(PYTHON) tools/make_exercises.py

check-exercises:
	$(PYTHON) tools/make_exercises.py --check

figures:
	$(PYTHON) tools/make_figures.py

format:
	clang-format -i --style="{BasedOnStyle: Google, IndentWidth: 4, ColumnLimit: 110}" \
		$$(git ls-files 'cpp/*.hpp' 'cpp/*.cpp' 'python/*.hpp' 'python/*.cpp')

clean:
	rm -rf build
