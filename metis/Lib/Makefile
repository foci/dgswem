SHELL:=/bin/sh

# CC:= pgcc

# CFLAGS:= -O2 -DLINUX 

OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o memory.o \
       pmetis.o pqueue.o refine.o util.o timing.o debug.o \
       bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
       kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
       mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
       myqsort.o compress.o parmetis.o estmem.o \
       mpmetis.o mcoarsen.o mmatch.o minitpart.o mbalance.o \
       mrefine.o mutil.o mfm.o mkmetis.o mkwayrefine.o mkwayfmh.o \
       mrefine2.o minitpart2.o mbalance2.o mfm2.o \
       kvmetis.o kwayvolrefine.o kwayvolfm.o subdomains.o

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

../libmetis.a: $(OBJS)
	ar -rv $@ $(OBJS)
	ar -ts $@
clean:
	rm -f *.o
