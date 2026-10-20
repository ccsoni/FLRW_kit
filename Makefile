CC       ?= gcc
CXX      ?= g++
CFLAGS   ?= -O2 -std=c11 -Wall -Wextra -pedantic
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic
AR       ?= ar
LDLIBS   += -lm

libflrw_kit.a: flrw_kit.o
	$(AR) rcs $@ $^

flrw_kit.o: flrw_kit.c flrw_kit.h
	$(CC) $(CFLAGS) -c -o $@ $<

test_flrw_kit: test_flrw_kit.c libflrw_kit.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lflrw_kit $(LDLIBS)

test_flrw_kit_cpp: test_flrw_kit_cpp.cpp flrw_kit.hpp flrw_kit.h libflrw_kit.a
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lflrw_kit $(LDLIBS)

test: test_c test_cpp

test_c: test_flrw_kit
	./test_flrw_kit

test_cpp: test_flrw_kit_cpp
	./test_flrw_kit_cpp

clean:
	rm -f flrw_kit.o libflrw_kit.a test_flrw_kit test_flrw_kit_cpp

.PHONY: test test_c test_cpp clean
