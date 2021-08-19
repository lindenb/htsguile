#ifndef BCFGUILE_H
#define BCFGUILE_H

#include <libguile.h>
#include "htslib/vcf.h"

#define COMMON_CTX_FIELDS 	char signature[4];\
	bcf_hdr_t *header;\
	bcf1_t *rec

typedef struct VariantContext {
	COMMON_CTX_FIELDS;
	}*VariantContextPtr;

typedef struct GenotypeContext {
	COMMON_CTX_FIELDS;
	int sample_idx;
	}GenotypeContext,*GenotypeContextPtr;

#define GENOTYPE_SIGNATURE "GTP"
#define VARIANT_SIGNATURE "VAR"

#undef COMMON_CTX_FIELDS

#endif

