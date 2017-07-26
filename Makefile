HEADERS := platypus/*.pxd
PYX := platypus/*.pyx
PY := platypus/*.py
FLAKE8=flake8 --max-line-length=120
SCRIPTS=bin/platypus

flake8:
	${FLAKE8} ${PY} test

cleanAll: clean cleanDocs
	rm -rf env

clean: cleanDocs
	pip uninstall -y Platypus
	find . -name __pycache__ | xargs rm -rf

cleanDocs:
	cd docs/sphinx; make clean

.PHONY:
pdfdocs: docs/sphinx/*.rst
	cd docs/sphinx; make latexpdf

.PHONY:
docs: docs/sphinx/*.rst
	cd docs/sphinx; make html;
	cp -rf docs/sphinx/_build/html/* docs/

wheels:
	pip wheel .

env/bin/platypus: ${HEADERS} ${PYX} ${PY} ${SCRIPTS}
	./install.sh
	touch env/bin/platypus

.PHONY:
install: env/bin/platypus ;

unittest: install
	@echo ''
	@echo 'Running unit tests'
	@echo ''
	@pytest test/unit

acceptancetest: install
	@echo ''
	@echo 'Running acceptance tests'
	@echo ''
	@pytest test/acceptance

smoketest: install
	@echo ''
	@echo 'Running smoke tests'
	@echo ''
	./test/smoke/check_installation_succeeded.bash

test: flake8 smoketest unittest acceptancetest
	@echo ''
	@echo 'Finished running all tests'
	@echo ''

test_coverage:
	pytest --cov=env/lib/python2.7/site-packages/platypus test
