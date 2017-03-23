PROGRAMS = strassen test
all: $(PROGRAMS)

include rules.mk

%.o: %.c $(BUILDSTAMP)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(DEPCFLAGS) $(O) -o $@ -c $<

strassen: strassen.o
	$(CC) $(CFLAGS) $(O) -o $@ $^ -lm

test: test.o
	$(CC) $(CFLAGS) $(O) -o $@ $^ -lm

clean:
	rm -f *.o $(PROGRAMS)
	rm -rf $(DEPSDIR) *.dSYM

.PHONY: all clean
