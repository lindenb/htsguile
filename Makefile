HTSLIB=../htslib
CC=gcc
CFLAGS=-g -Wall -I$(HTSLIB) `guile-config  compile`
LDFLAGS=-L$(HTSLIB) `guile-config  link`
bcfguile : bcfguile.o
	$(CC) -o $@   $^ $(LDFLAGS) -lhts

bcfguile.o : bcfguile.c bcfguile.h bcfguile.x
	$(CC) -c -o $@ $(CFLAGS) $<
bcfguile.x : bcfguile.c bcfguile.h
	guile-snarf -o $@ $(CFLAGS) $<

clean:
	rm -f bcfguile.x bcfguile.o bcfguile
