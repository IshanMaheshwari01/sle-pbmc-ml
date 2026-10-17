PY      ?= python
CONFIG  ?= config/config.yaml

.PHONY: help mock download inspect dev subsample qc integrate annotate phase1 phase1-dev test clean-processed

help:
	@echo "make mock       - build a small fake dataset and run the pipeline on it (no download)"
	@echo "make download   - download the 12 GB CELLxGENE h5ad (resumable)"
	@echo "make inspect    - data audit: design, confounding, counts check"
	@echo "make dev        - tiny subsample (40 cells/sample) + pseudobulk, for code iteration"
	@echo "make subsample  - full subsample (500 cells/sample) + pseudobulk over all cells"
	@echo "make qc         - phase 1: per-batch QC check on the subsample"
	@echo "make integrate  - phase 1: Harmony + scVI, benchmark, disease-signal and IFN checks"
	@echo "make annotate   - phase 1: CellTypist check of the author labels"
	@echo "make phase1     - all three phase-1 steps"
	@echo "make phase1-dev - all three on the tiny dev subsample (quick trial run)"
	@echo "make test       - run the test suite"

mock:
	$(PY) tests/mock_data.py data/raw/mock.h5ad
	$(PY) scripts/01_inspect.py --config config/mock.yaml
	$(PY) scripts/02_subsample.py --config config/mock.yaml

download:
	$(PY) scripts/00_download.py --config $(CONFIG)

inspect:
	$(PY) scripts/01_inspect.py --config $(CONFIG)

dev:
	$(PY) scripts/02_subsample.py --config $(CONFIG) --dev

subsample:
	$(PY) scripts/02_subsample.py --config $(CONFIG)

qc:
	$(PY) scripts/03_qc.py --config $(CONFIG)

integrate:
	$(PY) scripts/04_integrate.py --config $(CONFIG)

annotate:
	$(PY) scripts/05_annotate.py --config $(CONFIG)

phase1: qc integrate annotate

phase1-dev:
	$(PY) scripts/03_qc.py --config $(CONFIG) --dev
	$(PY) scripts/04_integrate.py --config $(CONFIG) --dev --max-epochs 20
	$(PY) scripts/05_annotate.py --config $(CONFIG) --dev

test:
	$(PY) -m pytest -q

clean-processed:
	rm -f data/processed/*.h5ad
