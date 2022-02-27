OBJ=pairlist.o
TESTS=tests/pairlist-test tests/pairlist-test2 tests/pairlist-test3
ALL=$(OBJ) pairlist-test2.data pairlist-test3.gro
PKGNAME=pairlist
GENICE=genice2

all: $(OBJ) README.md
	echo Done.

# Section: test

test: $(ALL)
	tests/pairlist-test
	tests/pairlist-test2
	tests/pairlist-test3
%.o: %.c
	$(CC) -c -g $< -o $@ -I.

%: %.c
%: %.o pairlist.o
	$(CC) -g $^ -o $@
pairlist-test2.data:
	$(GENICE) --format q --rep 3 3 3 1h --water tip4p > $@
pairlist-test3.gro: Makefile
	$(GENICE) --format g --rep 10 10 10  1h > $@


CRN1.gro:
	genice CRN1 -r 1 1 1 > $@
test2: CRN1.gro
	time python tests/test1_rdf.py < $< > test1.rdf #1.5 sec
	time python tests/test2_rdf.py < $< > test2.rdf #.26 sec
	time python tests/test3_rdf.py < $< > test3.rdf #.22 sec
CRN1x222.gro:
	genice CRN1 -r 2 2 2 > $@
test3: CRN1x222.gro
#	time python test1_rdf.py < $< > test1.rdf #95 sec
	time python tests/test20_rdf.py < $< > test20.rdf #3.6 sec
	time python tests/test2_rdf.py < $< > test2.rdf #1.2 sec
	time python tests/test3_rdf.py < $< > test3.rdf #0.8 sec

%: temp_% replacer.py pairlist.py benchmark/benchmark.py
	python3 replacer.py < $< > $@
	-fgrep '}}' $@


# Section: Deploy

test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install --no-cache-dir --index-url https://test.pypi.org/simple/ $(PKGNAME)


install:
	./setup.py install
uninstall:
	-pip uninstall -y pairlist
build: README.md
	./setup.py sdist # bdist_wheel


deploy: build
	twine upload --repository pypi dist/*
check:
	./setup.py check


# Section: Develop

pep8:
	autopep8 -r -a -a -i ./

clean:
	-rm $(ALL) *.so *~ */*~ *.o *.gro *.rdf
	-rm -rf build dist
	-rm -rf PairList.egg-info
