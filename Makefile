OBJ=pairlist.o tests/pairlist-test tests/pairlist-test2 tests/pairlist-test3
ALL=$(OBJ) pairlist-test2.data pairlist-test3.gro
PKGNAME=pairlist
all: $(ALL)
	tests/pairlist-test
	tests/pairlist-test2
	tests/pairlist-test3
%.o: %.c
	$(CC) -c -g $< -o $@ -I.
#for tests
%: %.c
%: %.o pairlist.o
	$(CC) -g $^ -o $@
%.rst: %.md
	md2rst $<
pairlist-test2.data:
	genice --format q --rep 3 3 3 1h --water tip4p > $@
pairlist-test3.gro: Makefile
	genice --format g --rep 10 10 10  1h > $@


test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PKGNAME)


install:
	./setup.py install
uninstall:
	-pip uninstall -y pairlist
build: README.md
	./setup.py sdist bdist_wheel


deploy: build
	twine upload --repository pypi dist/*
check:
	./setup.py check


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

clean:
	-rm $(ALL) *.so *~ */*~ *.o *.gro *.rdf
	-rm -rf build dist
	-rm -rf PairList.egg-info

