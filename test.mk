BAM?=test/toy.bam
FILTERSAM?=./htsguile filtersam
.PHONY:all

all: $(addprefix test,01 02 03)

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
	
