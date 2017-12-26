BAM?=test/toy.bam
FILTERSAM?=./htsguile filtersam
.PHONY:all

all: $(addprefix test,01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)

test01:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (and (> (hts-read-length R ) 10 ) (integer? (string-contains  (hts-read-seq R) "GATAAGGGATA")))))' $(BAM)
	
test02:
	echo '(define read-filter (lambda (R) (and (> (hts-read-length R ) 10 ) (integer? (string-contains  (hts-read-seq R) "GATAAGGGATA")))))' > $(addsuffix .scm,$@)
	$(FILTERSAM) -f   $(addsuffix .scm,$@) $(BAM)
	rm -f  $(addsuffix .scm,$@)
	
	
test03:
	$(FILTERSAM) -b -c 0 -e  '(define read-filter (lambda (R) (and (> (hts-read-length R ) 10 ) (integer? (string-contains  (hts-read-seq R) "GATAAGGGATA")))))' $(BAM) | file -

##  
test04:
	echo "(define read-filter (lambda (R) (= (hts-read-seq-at R 0) (hts-read-seq-at R 1))))" > $(addsuffix .scm,$@)
	$(FILTERSAM) -f   $(addsuffix .scm,$@) $(BAM)
	rm -f  $(addsuffix .scm,$@)
	
test05:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (= 83 (hts-read-flag R))))' $(BAM)
test06:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-paired? R)))' $(BAM)
test07:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-proper-pair? R)))' $(BAM)
test08:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-unmapped? R)))' $(BAM)
test09:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-mate-unmapped? R)))' $(BAM)	
test10:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-reverse-strand? R)))' $(BAM)	
test11:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-mate-reverse-strand? R)))' $(BAM)	
test12:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-1st-in-pair? R)))' $(BAM)	
test13:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-2nd-in-pair? R)))' $(BAM)	
test14:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-secondary? R)))' $(BAM)	
test15:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-qcfail? R)))' $(BAM)	
test16:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-duplicate? R)))' $(BAM)
test17:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-supplementary? R)))' $(BAM)
test18:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (> (hts-read-mapq R) 29)))' $(BAM)
test19:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (integer? (string-contains  (hts-read-cigar-string R) "I"))))' $(BAM)
test20:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-clipped? R)))' $(BAM)
test21:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (> (vector-length (hts-read-cigar R)) 2)))' $(BAM)
test22:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (string=? (hts-read-contig R) "ref")))' $(BAM)
test23:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (not (= (hts-read-tid R) 1))))' $(BAM)
test24:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (not (= (hts-read-pos R) (hts-read-unclipped-start R)))))' $(BAM)
test25:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (not (= (hts-read-end R) (hts-read-unclipped-end R)))))' $(BAM)
