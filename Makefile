HTSLIB=../htslib
CC=gcc
CFLAGS=-Wall -I$(HTSLIB)
LDFLAGS=-L$(HTSLIB)
bcfguile : bcfguile.o
	$(CC) -o $@  $(LDFLAGS) $^ -lhts 
bcfguile.o : bcfguile.c bcfguile.h
	$(CC) -c -o $@ $(CFLAGS) $<
