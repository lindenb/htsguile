BAM?=test/toy.bam
FILTERSAM?=./htsguile filtersam
.PHONY:all

all: $(addprefix test,01 02 03 04 05 06 07 08 09 10 11)

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
test11:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-1st-in-pair? R)))' $(BAM)	
test11:
	$(FILTERSAM) -e  '(define read-filter (lambda (R) (hts-read-2nd-in-pair? R)))' $(BAM)	
	
