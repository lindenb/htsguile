#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <libguile/snarf.h>
#include "bcfguile.h"

#define WHERE do { fprintf(stderr,"[%s][%d]",__FILE__,__LINE__);} while(0)
#define INFO(FormatLiteral, ...) do { WHERE;fprintf(stderr,"[INFO]" FormatLiteral "\n", ##__VA_ARGS__);} while(0)
#define WARNING(FormatLiteral, ...) do { WHERE;fprintf(stderr,"[WARNING]" FormatLiteral "\n", ##__VA_ARGS__);} while(0)
#define FATAL(FormatLiteral, ...) do { WHERE;fprintf(stderr,"[FATAL]" FormatLiteral "\n", ##__VA_ARGS__);abort()} while(0)
#define ASSERT_NOT_NULL(a) do {if((a)==NULL) { FATAL("null pointer") } while(0)

VariantContextPtr variant = NULL;

struct GenotypeShuttle {
	int phased;
	int nocall;
	int* alleles;
	size_t allele_count;
	size_t allele_capacity;
	int error_flag;
	int  max_ploidy;
	};


static void scanGenotype(VariantContextPtr ctx,int sample_idx, struct GenotypeShuttle* shuttle) {
	int32_t *gt_arr = NULL, ngt_arr = 0;
	bcf_unpack(ctx->rec,BCF_UN_IND);
	int j;
	int ngt = bcf_get_genotypes(ctx->header,ctx->rec, &gt_arr, &ngt_arr);
	int nsmpl = bcf_hdr_nsamples(ctx->header);
	int max_ploidy = ngt/nsmpl;
  	int32_t *ptr = gt_arr + max_ploidy*sample_idx;
  	shuttle->max_ploidy = max_ploidy;
  	for (j=0; j<max_ploidy; j++)
    	  {
         // if true, the sample has smaller ploidy
         if ( ptr[j]==bcf_int32_vector_end ) break;
  		// is phased?
  		if(bcf_gt_is_phased(ptr[j])) {
	         shuttle->phased = 1;
	         }
	         // the VCF 0-based allele index
	         int allele_index ;

	         // missing allele
	         if ( bcf_gt_is_missing(ptr[j]) ) {
	         	allele_index = -1;
	         	}
	         else
	         	{
		        allele_index= bcf_gt_allele(ptr[j]);
		        }
  			 if(shuttle->alleles != NULL) {
  			 	if(shuttle->allele_count + 1 >= shuttle->allele_capacity) {
  			 		shuttle->allele_capacity++;
  			 		shuttle->alleles = (int*)realloc(shuttle->alleles,shuttle->allele_capacity);
  			 		ASSERT_NOT_NULL(shuttle->alleles);
  			 		}
  			 	shuttle->alleles[shuttle->allele_count] = allele_index ;
  			 	}
  			shuttle->allele_count++;
  		    }
	}
/** dispose Genotype data created by scm_from_pointer in 'makeGenotype' */
static void freeGenotype(void* p) {
	GenotypeContextPtr ptr=(GenotypeContextPtr)p;
	if(ptr!=NULL) {
		free(ptr);
		}
	}

static int is_genotype(SCM g) {
	if(!scm_is_pointer(g)) return 0;
	GenotypeContextPtr ptr=(GenotypeContextPtr)scm_to_pointer(g);
	ASSERT_NOT_NULL(ptr);
	if(memcmp(ptr->signature,GENOTYPE_SIGNATURE,4*sizeof(char))!=0) return 0;
	return 1;
	}

/** create genotype from sample index */
static SCM makeGenotype(VariantContextPtr ctx,int sample_idx) {
	if(sample_idx<0 || sample_idx>=bcf_hdr_nsamples(ctx->header)) {
		FATAL("bad sample index 0<=%d<%d",sample_idx,bcf_hdr_nsamples(ctx->header));
		return SCM_UNDEFINED;
		}
	else
		{
		GenotypeContextPtr gt = (GenotypeContextPtr)malloc(sizeof(GenotypeContext));
		ASSERT_NOT_NULL(gt);
		memcpy((void*)(gt->signature),GENOTYPE_SIGNATURE,4*sizeof(char));
		gt->header = ctx->header;
		gt->rec = ctx->rec;
		gt->sample_idx = sample_idx;
		return scm_from_pointer((void*)gt,freeGenotype);
		}
	}

static int to_sample_index(VariantContextPtr ctx,SCM sexpgtidx) {
	int sample_idx=-1;
	if(is_genotype(sexpgtidx)) {
		sample_idx  = ((GenotypeContextPtr)scm_to_pointer)->sample_idx;
		}
	else if(scm_is_string(sexpgtidx)) {
		char* sn = scm_to_locale_string(sexpgtidx);
		if(sn!=NULL) {
			sample_idx = bcf_hdr_id2int(ctx->header,BCF_DT_SAMPLE,sn);
			if(sample_idx<0) {
				WHARNING("unknown sample \"%s\"",sn);
				}
			free(sn);
			}
		else
			{
			WHARNING("null sample");
			}
        }
	else if(scm_is_integer(sexpgtidx))
		{
		sample_idx = scm_to_int32(sexpgtidx);
		}
	else
		{
		FATAL("illegal argument: Sample is not a genotype, a string or an integer ");
		}
	return sample_idx;
	}

static GenotypeContextPtr* to_genotype_array(SCM param,int *nsamples) {
	GenotypeContextPtr *genotypes = NULL;
	if(scm_is_pair(param)) {
		int i=0;
		*nsamples= scm_to_int32(scm_length(param));
		genotypes=(GenotypeContextPtr*)malloc((*nsamples)*sizeof(GenotypeContextPtr));
		while(scm_is_pair(param)) {
			genotypes[i] =  SCM_CAR(param);
			param = SCM_CDR(param);
			i++;
			}
		return genotypes;
		}
	else if(is_genotype(param)) {
		genotypes=(GenotypeContextPtr*)malloc(sizeof(GenotypeContextPtr));
		genotypes[0] = makeGenotype(idx);
		return genotypes;
		}
	else
		{
		int idx = to_sample_index(param);
		*nsamples = 1;
		genotypes=(GenotypeContextPtr*)malloc(sizeof(GenotypeContextPtr));
		genotypes[0] = makeGenotype(idx);
		return genotypes;
		}
	}
	

static SCM VariantGetGenotype(VariantContextPtr ctx,SCM sexpgtidx) {
	return makeGenotype(to_sample_index(sexpgtidx));
	}

SCM_DEFINE (scm_bcf_has_variant_type, "is-variant-type?", 0, 1, 0,(SCM type),"test variant type") {
	char* s= scm_to_locale_string(type);
	int t=0;
	#define CMP_TYPE(T) if(strcasecmp(s,#T)==0) {t =  (bcf_get_variant_types(variant->rec) == VCF_##T);}
	CMP_TYPE(REF)
	else CMP_TYPE(SNP)
	else CMP_TYPE(MNP)
	else CMP_TYPE(INDEL)
	else CMP_TYPE(OTHER)
	else CMP_TYPE(BND)
	else CMP_TYPE(OVERLAP)
	free(s);
	#undef CMP_TYPE
	return scm_to_bool(t);
}

SCM_DEFINE (scm_bcf_variant_type, "variant-type", 0, 0, 1,(SCM ith_allele),"return variant type") {
int t;
if(SCM_IS_UNDETERMINDED(ith_allele)) {	
 t = bcf_get_variant_types(bcf1_t *rec);
 }
else if(scm_is_integer(ith_allele)) {
	t=bcf_get_variant_types(bcf1_t *rec,scm_to_int32(ith_allele));
	}
else
	{
	scm_wrong_type(ith_allele);
	}
#define CASE_T(X) case VCF_#X : return scm_from_locale_string(##X);break
switch(t)
	{
	CASE_T(REF)
	CASE_T(SNP)
	CASE_T(MNP)
	CASE_T(INDEL)
	CASE_T(OTHER)
	CASE_T(BND)
	CASE_T(OVERLAP)
	default: break;
	}
#undef CASE_T
return SCM_UNDEFINED;
}

SCM_DEFINE (scm_bcf_is_snp, "is-snp?", 0, 0, 0,(),"return true is variant is SNP") {
return scm_from_bool(bcf_is_snp(variant->rec));
}

SCM_DEFINE (scm_bcf_has_id, "has-id?", 0, 0, 0,(),"return true is variant has ID") {
return scm_from_bool(variant->rec->d.id);
}

SCM_DEFINE (scm_bcf_has_qual, "has-qual?", 0, 0, 0,(),"return true is variant has QUAL") {
return scm_from_bool(!bcf_float_is_missing(variant->rec->qual));
}

SCM_DEFINE (scm_bcf_chrom, "chrom", 0, 0, 0,(),"chromosome") {
return scm_from_locale_string(variant->header->id[BCF_DT_CTG][variant->rec->rid].key);
}

SCM_DEFINE (scm_bcf_pos, "pos", 0, 0, 0,  (), "pos") {
return scm_from_int64( variant->rec->pos + 1);
}

SCM_DEFINE (scm_bcf_id, "id", 0, 0, 0,(),"id") {
return scm_from_locale_string(variant->rec->d.id?variant->rec->d.id:"");
}



SCM_DEFINE (scm_bcf_nsamples, "n-samples", 0, 0, 0,  (), "count samples") {
return scm_from_int64( bcf_hdr_nsamples(variant->header));
}
SCM_DEFINE (scm_bcf_has_samples, "has-samples?", 0, 0, 0,  (), "has any sample") {
return scm_from_bool( bcf_hdr_nsamples(variant->header)>0);
}

SCM_DEFINE (scm_bcf_samples, "samples", 0, 0, 0,  (), "samples") {
int i;
SCM ret = SCM_UNDEFINED;
bcf_unpack(variant->rec,BCF_UN_ALL);
for(i=0; i<  bcf_hdr_nsamples(variant->header) ;i++) {
	SCM list = scm_list_1(scm_from_locale_string(variant->header->samples[i]));
	if(i==0) {
		ret = list;
		}
	else
		{
		ret = scm_append(scm_list_2(ret,list));
		}
	}
return ret;
}


SCM_DEFINE (scm_bcf_complement_samples, "not-samples", 1, 0, 0,  (SCM arg), "complement samples") {
	int i;
	SCM ret = SCM_UNDEFINED;
	bcf_unpack(variant->rec,BCF_UN_ALL);
	int* samples_indexes = NULL;
	int n_samples=0;

	for(i=0; i<  bcf_hdr_nsamples(variant->header) ;i++) {
		if(scm_is_pair(arg)) {
			int sample_index = -1;
			SCM list = arg;
			while(scm_is_pair(param)) {
				SCM car =  SCM_CAR(list);
				sample_index = to_sample_index(variant,car);
				if (sample_index == i) break;
				list = SCM_CDR(list);
				}
			if (sample_index == i) continue; 
			}
		else
			{
			int idx = to_sample_index(param);
			if (gt->sample_index == i) continue; 
			}


	SCM list = scm_list_1(scm_from_locale_string(variant->header->samples[i]));
	if(scm_is_undefined(ret)) {
		ret = list;
		}
	else
		{
		ret = scm_append(scm_list_2(ret,list));
		}
	}
return ret;
}

SCM_DEFINE (scm_bcf_nalleles, "n-alleles", 0, 0, 0,  (), "count alleles") {
return scm_from_int64( variant->rec->n_allele);
}
SCM_DEFINE (scm_bcf_allele, "allele", 1, 0, 0,  (SCM arg_idx), "allele") {
bcf_unpack(variant->rec,BCF_UN_ALL);
if(!scm_is_integer(arg_idx)) {
	scm_wrong_type_arg("not an integer",1,arg_idx);
	}
int idx = scm_to_int32(arg_idx);
if(idx<0 || idx>=variant->rec->n_allele) scm_wrong_type_arg("out of range",1,arg_idx);
return scm_from_locale_string(variant->rec->d.allele[idx]);
}

SCM_DEFINE (scm_bcf_alleles, "alleles", 0, 0, 0,  (), "alleles") {
int i;
SCM ret = SCM_UNDEFINED;
bcf_unpack(variant->rec,BCF_UN_ALL);
for(i=0; i<  variant->rec->n_allele;i++) {
	SCM list = scm_list_1(scm_from_locale_string(variant->rec->d.allele[i]));
	if(i==0) {
		ret = list;
		}
	else
		{
		ret = scm_append(scm_list_2(ret,list));
		}
	}
return ret;
}

SCM_DEFINE (scm_bcf_qual, "qual", 0, 0, 0,  (), "qual") {
if(bcf_float_is_missing(variant->rec->qual)) return SCM_UNDEFINED;
return scm_from_double( variant->rec->qual );
}

SCM_DEFINE (scm_bcf_has_filter, "has-filter?", 1, 0, 0,(SCM arg_filter),"return true is variant has filter") {
if(scm_is_pair(arg_filter)) {
	while(scm_is_pair(arg_filter)) {
		SCM p = CAR(arg_filter);
		if(!scm_is_string(p))  scm_wrong_type_arg("not a string",1,p);
		char *s=scm_to_locale_string(p);
		int f= bcf_has_filter(variant->header,variant->rec,s);
		free(s);
		if(f) return SCM_TRUE;
		arg_filter = CDS(arg_filter);
		}
	return SCM_FALSE;
	}
else
	{
	if(!scm_is_string(arg_filter))  scm_wrong_type_arg("not a string",1,arg_filter);
	char *s=scm_to_locale_string(arg_filter);
	SCM ret = scm_from_bool(!bcf_has_filter(variant->header,variant->rec,s));
	free(s);
	return ret;
	}
}

SCM_DEFINE (scm_bcf_is_filtered, "filtered?", 0, 0, 0,  (), "is filtered") {
int i;
SCM ret = SCM_UNDEFINED;
bcf_unpack(variant->rec,BCF_UN_ALL);
if( variant->rec->d.n_flt==0) return SCM_FALSE;
// doc: "PASS" and "." can be used interchangeably
if(bcf_has_filter(variant->header,variant->rec,"PASS")) return SCM_FALSE;
return SCM_TRUE;
}

SCM_DEFINE (scm_bcf_filters, "filters", 0, 0, 0,  (), "filters") {
int i;
SCM ret = SCM_UNDEFINED;
bcf_unpack(variant->rec,BCF_UN_ALL);
for(i=0; i<  variant->rec->d.n_flt;i++) {
	SCM list = scm_list_1(scm_from_locale_string(variant->header->id[BCF_DT_ID][variant->rec->d.flt[i]].key));
	if(i==0) {
		ret = list;
		}
	else
		{
		ret = scm_append(scm_list_2(ret,list));
		}
	}
return ret;
}

SCM_DEFINE (scm_bcf_genotypes, "genotypes", 0, 0, 0,  (), "genotypes") {
int i;
SCM ret = SCM_UNDEFINED;
bcf_unpack(variant->rec,BCF_UN_ALL);
for(i=0; i<  bcf_hdr_nsamples(variant->header) ;i++) {
	SCM list = scm_list_1(makeGenotype(i));
	if(i==0) {
		ret = list;
		}
	else
		{
		ret = scm_append(scm_list_2(ret,list));
		}
	}
return ret;
}

SCM_DEFINE (scm_bcf_genotype, "genotype", 1, 0, 0,  (SCM sampleoridx), "genotype") {
if(scm_is_pair(sampleoridx)) {
	int i=0;
	SCM ret= SCM_UNDEFINED;
	while(scm_is_pair(sampleoridx)) {
		SCM x = SCM_CAR(sampleoridx);
		SCM list = scm_list_1(VariantGetGenotype(variant,x));
		if( i == 0) {
			ret = list;
		} else {
			ret = scm_append(scm_list_2(ret,list));
		}
		sampleoridx = SCM_CDR(sampleoridx);
		i++;
		}
	return ret;
	}
else
	{
	return VariantGetGenotype(variant,sampleoridx);
	}
}

SCM_DEFINE (scm_genotype_is_phased, "gt-is-phased?", 0, 0, 0,  (SCM genotype), "genotype is phased") {
	scanGenotype(VariantContextPtr ctx,int sample_idx, struct GenotypeShuttle* shuttle)
	return SCM_UNDEFINED;
}
SCM_DEFINE (scm_genotype_is_hom, "gt-is-hom?", 0, 0, 0,  (SCM genotype), "genotype is homozygous") {
	SCM ret = SCM_UNDEFINED;
	if(!scm_is_pair(genotype)) {
		return  scm_genotype_is_hom(scm_list_1(genotype));
		}
	while(scm_is_pair(genotype)) {
		SCM p = CAR(genotype);
		CALL(fun,p,val);
		SCM list = scm_list_1(val);
		if(scm_is_undefined()) {
			SCM val;
			ret = list;
			}
		else
			{
			scm_append(scm_list_2(ret,list));
			}
		genotype = CDS(genotype);
		}
	return ret;
	}
	
SCM_DEFINE (scm_genotype_is_hom_ref, "gt-is-hom-ref?", 0, 0, 0,  (SCM genotype), "genotype is homozygous") {
	return SCM_UNDEFINED;
}
SCM_DEFINE (scm_genotype_is_hom_var, "gt-is-hom-var?", 0, 0, 0,  (SCM genotype), "genotype is homozygous") {
	return SCM_UNDEFINED;
}
SCM_DEFINE (scm_genotype_is_het, "gt-is-het?", 0, 0, 0,  (SCM genotype), "genotype is homozygous") {
	return SCM_UNDEFINED;
}
SCM_DEFINE (scm_genotype_alleles, "gt-alleles", 1, 0, 0,  (SCM param), "genotype alleles") {
	GenotypeShuttle shuttle;
	SCM ret = SCM_UNDEFINED;
	memset((void*)&shuttle,0,sizeof(GenotypeShuttle));
	shuttle.allele_capacity=10;
	shuttle.alleles = (int*)malloc(shuttle.allele_capacity*sizeof(int));
	int sample_index= to_sample_index(param);
	scanGenotype(variant,sample_index,&shuttle);
	if(shuttle.error_flag)
		{
		ret = SCM_UNDEFINED;
		}
	else
		{
		int i=0;
		for(i=0;i< shuttle.allele_count;i++) {
		SCM list = scm_list_1(scm_from_int32(shuttle.alleles[i]));
		if( i == 0) {
			ret = list;
		} else {
			ret = scm_append(scm_list_2(ret,list));
			}
		}
		}
	free(shuttle.alleles);
	return ret;
}


#ifdef xxx
static SCM fmt_array(int n, int type, void *data)
{
    int j = 0;
    uint32_t e = 0;
    
    if (type == BCF_BT_CHAR)
    {
    	SCM list = SCM_EOL;
    	kstring_t ks;
        char *p = (char*)data;
        for (j = 0; j < n && *p; ++j, ++p)
        	{
            if ( *p==bcf_str_missing ) {
            	list = scm_append(scm_list_2(SCM_UNDEFINED,list));
            	}
            else if (*p==',') {
            	SCM lcls = scm_from_locale_string(ks.s);
            	list = scm_append(scm_list_2(list,lcls));
            	}
            else {
            	kputc(*p, s);
            	}
        	}
        SCM lcls = scm_from_locale_string(ks.s);
        list = scm_append(scm_list_2(lcls,list));
        return list;
    }
    else
    {
    	SCM list = SCM_EOL;
        #define BRANCH(type_t, convert, is_missing, is_vector_end, convert2) { \
            uint8_t *p = (uint8_t *) data; \
            for (j=0; j<n; j++, p += sizeof(type_t))    \
            { \
                type_t v = convert(p); \
                if ( v == is_vector_end ) break; \
                if ( v == is_missing )  list = scm_append(scm_list_2(SCM_UNDEFINED,list)); \
                else {\
                	SCM value = convert ;\
                	list = scm_append(scm_list_2(value,list)); \
                	}\
            } \
            return list; \
        }
        switch (type) {
            case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8, bcf_int8_missing, bcf_int8_vector_end, scm_from_int8); break;
            case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, bcf_int16_missing,bcf_int16_vector_end, scm_from_int16); break;
            case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, bcf_int32_missing,bcf_int32_vector_end, scm_from_int32); break;
            case BCF_BT_FLOAT: BRANCH(uint32_t, le_to_u32, bcf_float_missing,bcf_float_vector_end, scm_from_double(le_to_float(p),s)); break;
            default: hts_log_error("Unexpected type %d", type); exit(1); break;
        }
        #undef BRANCH
    }
   return list;
}
#endif

SCM_DEFINE (scm_bcf_attributes, "attributes", 0, 0, 0,  (), "INFO attributes") {
int i;
SCM alist = SCM_EOL;
bcf_unpack(variant->rec,BCF_UN_ALL);
for (i = 0; i < variant->rec->n_info; ++i) {
    bcf_info_t *z = &variant->rec->d.info[i];
    if ( !z->vptr ) continue;
   
    SCM scm_key = scm_from_locale_string(variant->header->id[BCF_DT_ID][z->key].key);
    SCM scm_value = SCM_UNDEFINED;
    if (z->len <= 0) {
    	scm_value = SCM_BOOL_T;
    	}
    else
    	{
    	scm_value = SCM_UNDEFINED;//TODO fix this
    	}
    alist = scm_acons(scm_key,scm_value,alist);
	}
return alist;
}

int main(int argc,char** argv) {
int c;
char* filenameout = NULL;
char* guilexpr = NULL;
struct VariantContext ctx;
variant = &ctx;
while ((c = getopt (argc, argv, "o:e:")) != -1)
 {
 switch (c)
    {
    case 'e': guilexpr = optarg; break;
    case 'o': filenameout = optarg; break;
	case '?':
       	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
       	return EXIT_FAILURE;
    default:
    	fprintf (stderr, "Bad input.\n");
       	return EXIT_FAILURE;
    }
 }
/* initialize guile */
scm_init_guile();
#ifndef SCM_MAGIC_SNARFER
#include "bcfguile.x"
#endif
 
 
if(guilexpr==NULL) {
	fprintf(stderr,"Undefined scheme expression.\n");
	return EXIT_FAILURE;
	}
	

	
if(!(optind==argc || optind+1==argc)) {
	fprintf(stderr,"Illegal number of arguments.\n");
	return EXIT_FAILURE;
	}
htsFile* vcfin = hts_open(optind==argc?"-":argv[optind],"r");
if(vcfin==NULL) {
	fprintf(stderr,"Cannot open \"%s\". %s\n",optind==argc?"<stdin>":argv[optind],strerror(errno));
	return EXIT_FAILURE;
	}
ctx.header= bcf_hdr_read(vcfin);
if(ctx.header==NULL) {
	fprintf(stderr,"Cannot read header from \"%s\".\n",optind==argc?"<stdin>":argv[optind]);
	return EXIT_FAILURE;
	}
ctx.rec = bcf_init1();

htsFile* vcfout = hts_open(filenameout==NULL?"-":filenameout,"w");
if(vcfout==NULL) {
	fprintf(stderr,"Cannot write \"%s\". %s\n",filenameout==NULL?"<stdout>":filenameout,strerror(errno));
	return EXIT_FAILURE;
	}
if(bcf_hdr_write(vcfout,ctx.header)!=0) {
	fprintf(stderr,"Cannot write header.\n");
	return EXIT_FAILURE;
	}
while ( bcf_read(vcfin, ctx.header, ctx.rec)==0 ) {
	SCM ret = scm_c_eval_string (guilexpr);
	if ( !scm_is_bool(ret)) {
		fprintf(stderr,"Returned guile value is not a boolean.\n");
		return EXIT_FAILURE;
		}
	if(scm_is_true(ret)) {
		if ( bcf_write(vcfout,ctx.header,ctx.rec)!=0 ) {
			fprintf(stderr,"I/O error cannot write record.\n");
			return EXIT_FAILURE;
			}
		}
	}
bcf_hdr_destroy(ctx.header);
bcf_destroy(ctx.rec);
hts_close(vcfin);
hts_close(vcfout);
return EXIT_SUCCESS;
}
