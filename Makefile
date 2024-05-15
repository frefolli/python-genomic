@all:
	python -m build --wheel --no-isolation

install:
	python -m installer --destdir="${DESTDIR}" dist/*.whl
	mkdir -p ${DESTDIR}/usr/share/python-genomic
	cp -r samples ${DESTDIR}/usr/share/python-genomic

run:
	./actions.sh run

sonarqube:
	./actions.sh sonarqube

test:
	./actions.sh test

CI:
	./actions.sh CI

verify:
	./actions.sh verify
