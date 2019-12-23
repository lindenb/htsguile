%.o:%.c
	gcc  -o $@ -c  `guile-config compile` -I ${HOME}/src/htslib $<

test: ./a.out test.scm
	 LD_LIBRARY_PATH=/home/lindenb/src/htslib ./a.out < test.scm

./a.out : scm_htsfile.o scm_bcf_hdr.o scm_bcf_rec.o test.o
	gcc  -o $@ -L ${HOME}/src/htslib  $^  `guile-config link` -lhts

clean:
