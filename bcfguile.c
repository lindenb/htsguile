#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <libguile/snarf.h>
#include "bcfguile.h"

#define WHERE do { fprintf(stderr,"[%s][%d]\n",__FILE__,__LINE__);} while(0)

VariantContextPtr variant = NULL;

SCM_DEFINE (scm_bcf_is_snp, "is-snp?", 0, 0, 0,(),"return true is variant is SNP") {
return scm_from_bool(bcf_is_snp(variant->rec));
}

SCM_DEFINE (scm_bcf_has_id, "has-id?", 0, 0, 0,(),"return true is variant has ID") {
return scm_from_bool(variant->rec->d.id);
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
