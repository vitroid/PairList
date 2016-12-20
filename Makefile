OBJ=pairlist.o pairlist-test pairlist-test2 pairlist-test3
ALL=$(OBJ) pairlist-test2.data pairlist-test3.gro
all: $(ALL)
%.o: %.c
	$(CC) -c -g $< -o $@
#for tests
%: %.c
%: %.o pairlist.o
	$(CC) -g $^ -o $@
pairlist-test2.data:
	( cd ../GenIce; ./genice --format q --rep 3 3 3 1h --water tip4p ) > $@
pairlist-test3.gro:
	( cd ../GenIce; ./genice --format g --rep 3 3 3 1h ) > $@

clean:
	-rm $(OBJ)
