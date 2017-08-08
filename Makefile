HEADERS := platypus/*.pxd
PYX := platypus/*.pyx
PY := platypus/*.py
FLAKE8=flake8 --max-line-length=120
SCRIPTS=bin/platypus

cleanAll: clean
	rm -rf env

clean:
	pip uninstall -y Platypus
	find . -name __pycache__ | xargs rm -rf

.PHONY:
pdfdocs: docs/sphinx/*.rst
	cd docs/sphinx; make latexpdf

.PHONY:
docs: docs/sphinx/*.rst
	cd docs/sphinx; make html;
	cp -rf docs/sphinx/_build/html/* docs/
	cd docs/sphinx; make clean

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

flake8:
	${FLAKE8} ${PY} test
