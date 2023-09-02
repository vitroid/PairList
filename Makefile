OBJ=pairlist.o
TESTS=tests/pairlist-test tests/pairlist-test2 tests/pairlist-test3
ALL=$(OBJ)
PKGNAME=pairlist
GENICE=genice2

all: $(OBJ) README.md
	echo Done.

# Section: test

test: $(ALL) test1
	echo Done.
	# tests/pairlist-test
	# tests/pairlist-test2
	# tests/pairlist-test3
%.o: %.c
	$(CC) -c -g $< -o $@ -I.

%: %.c
%: %.o pairlist.o
	$(CC) -g $^ -o $@
# pairlist-test2.data:
# 	$(GENICE) --format q --rep 3 3 3 1h --water tip4p > $@
# pairlist-test3.gro: Makefile
# 	$(GENICE) --format g --rep 10 10 10  1h > $@


CRN1x222.gro:
	genice2 CRN1 -r 2 2 2 > $@
test1: CRN1x222.gro
	time python samples/HB.py < $<
	time python samples/RDF_crude.py < $< > RDF_crude.txt
	time python samples/RDF_pairlist.py < $< > RDF_pairlist.txt

%: temp_% replacer.py pairlist.py benchmark/benchmark.py
	python3 replacer.py < $< > $@
	-fgrep '}}' $@


# Section: Deploy

requirements.txt:
	pipenv lock -r > $@

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
