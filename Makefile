CC      ?= gcc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -pedantic
AR      ?= ar
LDLIBS  += -lm

libflrw_kit.a: flrw_kit.o
	$(AR) rcs $@ $^

flrw_kit.o: flrw_kit.c flrw_kit.h
	$(CC) $(CFLAGS) -c -o $@ $<

test_flrw_kit: test_flrw_kit.c libflrw_kit.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lflrw_kit $(LDLIBS)

test: test_flrw_kit
	./test_flrw_kit

clean:
	rm -f flrw_kit.o libflrw_kit.a test_flrw_kit

.PHONY: test clean
