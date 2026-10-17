PY      ?= python
CONFIG  ?= config/config.yaml

.PHONY: help mock download inspect dev subsample test clean-processed

help:
	@echo "make mock       - build a small fake dataset and run the pipeline on it (no download)"
	@echo "make download   - download the 12 GB CELLxGENE h5ad (resumable)"
	@echo "make inspect    - data audit: design, confounding, counts check"
	@echo "make dev        - tiny subsample (40 cells/sample) + pseudobulk, for code iteration"
	@echo "make subsample  - full subsample (500 cells/sample) + pseudobulk over all cells"
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

test:
	$(PY) -m pytest -q

clean-processed:
	rm -f data/processed/*.h5ad
