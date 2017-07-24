
PYTHON := python2.6
HEADERS := src/cython/cgenotype.pxd src/cython/cfilter.pxd
SOURCES := src/cython/cpopulation.pyx src/cython/cgenotype.pyx src/cython/cwindow.pyx src/cython/cfilter.pyx

all: ${HEADERS} ${SOURCES} cortex
	echo 'Building Platypus'
	echo 'C Libraries will be installed to ' ${PLATYPUS}
	cd src/cython; ${PYTHON} setup.py install --prefix=${PLATYPUS}
	cd src/python; ${PYTHON} setup.py install --prefix=${PLATYPUS}

platypus: ${HEADERS} ${SOURCES}
	echo 'Building Platypus'
	echo 'C Libraries will be installed to ' ${PLATYPUS}
	cd src/cython; ${PYTHON} setup.py install --prefix=${PLATYPUS}
	cd src/python; ${PYTHON} setup.py install --prefix=${PLATYPUS}

cortex:
	echo 'Building Cortex'
	cd src/cortex_var; make cortex_var 32_BITS=1

clean:
	cd src/cortex_var; make clean
	cd src/cython; rm -rf build; rm cpopulation.c cgenotype.c cwindow.c cfilter.c
	echo ${PLATYPUS}
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cpopulation.so
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cgenotype.so
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cwindow.so
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cfilter.so
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/variantcaller.so
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/variantFilter.so
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/platypusexceptions.py
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/window.py
	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/variantutils.py
