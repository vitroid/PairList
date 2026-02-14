OBJ=pairlist.o
TESTS=tests/pairlist-test tests/pairlist-test2 tests/pairlist-test3
ALL=$(OBJ)
PKGNAME=pairlist
GENICE=genice2

# Standalone C library for use from Nim, Julia, C++, etc.
CSOURCE=csource
LIBPAIRLIST=libpairlist.a

all: $(OBJ) README.md
	echo Done.

# Build static C library (libpairlist.a)
lib: $(LIBPAIRLIST)
$(LIBPAIRLIST): $(CSOURCE)/pairlist.o
	$(AR) rcs $@ $^
$(CSOURCE)/pairlist.o: $(CSOURCE)/pairlist.c $(CSOURCE)/pairlist.h
	$(CC) -c $(CFLAGS) -std=c99 -I$(CSOURCE) -o $@ $<

# Section: test

test: $(ALL) test1
	echo Done.
%.o: %.c
	$(CC) -c -g $< -o $@ -I.

%: %.c
%: %.o pairlist.o
	$(CC) -g $^ -o $@


CRN1x222.gro:
	genice2 CRN1 -r 2 2 2 > $@
test1: CRN1x222.gro
	time python samples/HB.py < $<
	time python samples/RDF_crude.py < $< > RDF_crude.txt
	time python samples/RDF_pairlist.py < $< > RDF_pairlist.txt

# https://qiita.com/yukinarit/items/0996180032c077443efb
# https://zenn.dev/atu4403/articles/python-githubpages
doc: README.md ./pairlist.py # CITATION.cff 
	pdoc -o docs ./pairlist.py --docformat google


%: temp_% replacer.py pairlist.py benchmark.ipynb pyproject.toml
	python3 replacer.py < $< > $@


# Section: Deploy

requirements.txt:
	pipenv lock -r > $@

test-deploy:
	poetry publish --build -r testpypi
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PKGNAME)
uninstall:
	-pip uninstall -y $(PKGNAME)
build: README.md $(wildcard cycles/*.py)
	poetry build
deploy:
	poetry publish --build
check:
	poetry check



# Section: Develop

pep8:
	autopep8 -r -a -a -i ./

clean:
	-rm $(ALL) *.so *~ */*~ *.o *.gro *.rdf
	-rm -f $(LIBPAIRLIST) $(CSOURCE)/pairlist.o
	-rm -rf build dist
	-rm -rf PairList.egg-info
	-rm setup.py # Generated automatically by Poetry
