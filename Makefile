OBJ=pairlist.o pairlist-test pairlist-test2 pairlist-test3
ALL=$(OBJ) pairlist-test2.data pairlist-test3.gro
all: $(ALL)
	./pairlist-test
	./pairlist-test2
	./pairlist-test3
#	python3 example.py
%.o: %.c
	$(CC) -c -g $< -o $@
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
register:
	./setup.py register -r pairlist
pypi:
	make README.rst
	./setup.py check
	./setup.py sdist bdist_wheel upload
install:
	make README.rst
	./setup.py install
build.:
	-rm *.so
	-rm -rf build
	python setup.py build_ext --inplace

CRN1.gro:
	genice CRN1 -r 1 1 1 > $@
test2: CRN1.gro
	time python test1_rdf.py < $< > test1.rdf #1.5 sec
	time python test2_rdf.py < $< > test2.rdf #.26 sec
	time python test3_rdf.py < $< > test3.rdf #.22 sec
CRN1x222.gro:
	genice CRN1 -r 2 2 2 > $@
test3: CRN1x222.gro
#	time python test1_rdf.py < $< > test1.rdf #95 sec
	time python test2_rdf.py < $< > test2.rdf #1.2 sec
	time python test3_rdf.py < $< > test3.rdf #0.8 sec

clean:
	-rm $(ALL) *.so *~ */*~ *.o *.gro *.rdf
	-rm -rf build dist
	-rm -rf PairList.egg-info

