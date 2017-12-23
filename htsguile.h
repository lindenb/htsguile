#ifndef HTSGUILE_H
#define HTSGUILE_H

typedef struct hts_guile_context_t {
	const bam_hdr_t *header;
	bam1_t *b;
	} HtsGuileCtx,*HtsGuileCtxPtr; 
	
	
#endif

