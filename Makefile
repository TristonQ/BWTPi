# Makefile produced by _mmakefile

CC = gcc
LDFLAGS =
CFLAGS = -std=gnu99 -Wall -ggdb

bencode : bencode.o
	$(CC) $(CFLAGS) -o bencode bencode.o $(LDFLAGS)
	
pencode : pencode.o
	$(CC) $(CFLAGS) -o pencode pencode.o $(LDFLAGS)
	
pencode2 : pencode2.o
	$(CC) $(CFLAGS) -o pencode2 pencode2.o $(LDFLAGS)
	
psearch : psearch.o
	$(CC) $(CFLAGS) -o psearch psearch.o $(LDFLAGS)

