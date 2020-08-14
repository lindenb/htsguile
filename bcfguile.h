#ifndef BCFGUILE_H
#define BCFGUILE_H

#include "htslib/vcf.h"

typedef struct VariantContext {
	bcf_hdr_t *header;
	bcf1_t *rec;
	}*VariantContextPtr;

#endif

