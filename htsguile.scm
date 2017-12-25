(define hts-read-clipped? (lambda (R) (or 
  (integer? (string-contains  (hts-read-cigar-string R) "S"))
  (integer? (string-contains  (hts-read-cigar-string R) "H"))
  )))
