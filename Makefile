OBJ=pairlist.o pairlist-test pairlist-test2 pairlist-test3
ALL=$(OBJ) pairlist-test2.data pairlist-test3.gro
all: $(ALL)
	./pairlist-test
	./pairlist-test2
	./pairlist-test3
#	python3 example.py
%.o: %.c
	$(CC) -c -g -O $< -o $@
#for tests
%: %.c
%: %.o pairlist.o bst.o
	$(CC) -g -O $^ -o $@
pairlist-test2.data:
	genice --format q --rep 3 3 3 1h --water tip4p > $@
pairlist-test3.gro: Makefile
	genice --format g --rep 20 20 20  1h > $@
TPPI.ar3r:
	genice TPPI --rep 1 1 1 --format r > TPPI.ar3r
6.ar3r:
	genice 6 --rep 1 1 1 --format r > 6.ar3r
7.ar3r:
	genice 7 --rep 1 1 1 --format r > 7.ar3r
%.TPPI.match: %.gro TPPI.ar3r
	./matcher $< TPPI.ar3r 0.02 0.25 > $@
%.TPPI.match.yap: %.TPPI.match TPPI.ar3r
	python3 ./match2yap.py $*.gro TPPI.ar3r < $< > $@
%.6.match: %.gro 6.ar3r
	./matcher $< 6.ar3r 0.02 0.25 > $@
%.6.match.yap: %.6.match 6.ar3r
	python3 ./match2yap.py $*.gro 6.ar3r < $< > $@
%.7.match: %.gro 7.ar3r
	./matcher $< 7.ar3r 0.02 0.25 > $@
%.7.match.yap: %.7.match 7.ar3r
	python3 ./match2yap.py $*.gro 7.ar3r < $< > $@
clean:
	-rm $(ALL)
