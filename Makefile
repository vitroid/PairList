all: pairlist.o pairlist-test pairlist-test.data
%.o: %.c
	$(CC) -c -g $< -o $@ -lm
pairlist-test.data:
	( cd ../GenIce; ./genice --format q --rep 3 3 3 1h --water tip4p ) > $@
pairlist-test: pairlist-test.o pairlist.o
	$(CC) -g $^ -o $@ -lm

