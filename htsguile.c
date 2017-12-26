#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <assert.h>

#include <libguile.h>
#include <libguile/modules.h> 
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htsguile.h"
#include "htsguile.scm.h"

#define UNUSED
#define GUILE_FILTER_NAME "read-filter" 

typedef struct cigaroperator {
  char name;
  int consume_read;
  int consume_ref;
  } CigarOperator,*CigarOperatorPtr;

static CigarOperator CIGAR_OP_M={'M',1,1};
static CigarOperator CIGAR_OP_I={'I',1,0};
static CigarOperator CIGAR_OP_S={'S',1,0};
static CigarOperator CIGAR_OP_H={'H',0,0};
static CigarOperator CIGAR_OP_D={'D',0,1};
static CigarOperator CIGAR_OP_N={'N',0,1};
static CigarOperator CIGAR_OP_X={'X',1,1};
static CigarOperator CIGAR_OP_EQ={'=',1,1};
static CigarOperator CIGAR_OP_P={'P',0,0};

static CigarOperatorPtr CigarOperatorFromName(char c)
  {
  switch(c)
    {
    case 'M': return &CIGAR_OP_M;
    case 'I': return &CIGAR_OP_I;
    case 'S': return &CIGAR_OP_S;  
    case 'H': return &CIGAR_OP_H;
    case 'D': return &CIGAR_OP_D;  
    case 'N': return &CIGAR_OP_N;
    case 'X': return &CIGAR_OP_X;
    case '=': return &CIGAR_OP_EQ;
    case 'P': return &CIGAR_OP_P;
    default: fprintf(stderr,"BAD CIGAR OPERATOR %c\n",c);exit(EXIT_FAILURE);
    }
  }

typedef struct cigarelement {
  int length;
  CigarOperatorPtr op;
  } CigarElement,*CigarElementPtr;


static CigarElementPtr makeCigarElement(int length,char op)
  {
  CigarElementPtr ce = (CigarElementPtr)malloc(sizeof(CigarElement));
  ce->length = length;
  ce->op = CigarOperatorFromName(op);
  return ce;
  }

static SCM make_mod = SCM_EOL;

static HtsGuileCtxPtr cast_to_ctx_ptr(SCM scm_ctx)
  {
  HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
  return ptr;
  }

static SCM hts_read_name(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return  scm_from_locale_string(bam_get_qname(ptr->b));
	}

static SCM hts_read_tid(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.tid);
	}

static SCM hts_read_mapq(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return  scm_from_signed_integer((int)ptr->b->core.qual);
	}

static SCM hts_read_contig(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
    	if ( ptr->header !=NULL &&
    	     ptr->b->core.tid >=0 &&
    	     ptr->b->core.tid < ptr->header->n_targets
    	     )
    		{
    		scm_from_locale_string(ptr->header->target_name[ptr->b->core.tid]);
    		}
    	return SCM_UNDEFINED;
    	}

static SCM hts_read_pos(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.pos + 1);
	}
	
static SCM hts_read_unclipped_start(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	if(ptr->b->core.flag & BAM_FUNMAP) return scm_from_signed_integer(0);
	int pos = ptr->b->core.pos + 1;
  int i;
	uint32_t *cigar = bam_get_cigar(ptr->b);
  for (i = 0; i < ptr->b->core.n_cigar; ++i) 
      {
      switch(bam_cigar_opchr(cigar[i]))
        {
        case 'S':
        case 'H': pos -= bam_cigar_oplen(cigar[i]);break;
        default: return scm_from_signed_integer(pos);
        }
      }
  return  scm_from_signed_integer(pos);
	}

static SCM hts_read_reference_length(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	if(ptr->b->core.flag & BAM_FUNMAP) return scm_from_signed_integer(0);
  int i,L=0;
	uint32_t *cigar = bam_get_cigar(ptr->b);
  for (i = 0; i < ptr->b->core.n_cigar; ++i) 
      {
      switch(bam_cigar_opchr(cigar[i]))
        {
        case 'M':
        case 'D':
        case 'N':
        case 'X':
        case '=': L+= bam_cigar_oplen(cigar[i]);break;
        default: break;
        }
      }
  return  scm_from_signed_integer(L);
	}

static SCM hts_read_reference_end(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	if(ptr->b->core.flag & BAM_FUNMAP) return scm_from_signed_integer(0);
  int i=0,pos = ptr->b->core.pos + 1 ;
	uint32_t *cigar = bam_get_cigar(ptr->b);
  for (i = 0; i < ptr->b->core.n_cigar; ++i) 
      {
      switch(bam_cigar_opchr(cigar[i]))
        {
        case 'M':
        case 'D':
        case 'N':
        case 'X':
        case '=': pos+= bam_cigar_oplen(cigar[i]);break;
        default: break;
        }
      }
  return  scm_from_signed_integer(pos);
	}

static SCM hts_read_unclipped_end(SCM scm_ctx)
	{
  HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	if(ptr->b->core.flag & BAM_FUNMAP) return scm_from_signed_integer(0);
  int pos = ptr->b->core.pos + 1 ;
	uint32_t *cigar = bam_get_cigar(ptr->b);
	int i=0;
	while(i < ptr->b->core.n_cigar) 
      {
      char op=bam_cigar_opchr(cigar[i]);
      if(!(op=='S' || op=='H')) break;
      i++;
      }
  while(i < ptr->b->core.n_cigar) 
      {
      switch(bam_cigar_opchr(cigar[i]))
        {
        case 'H':
        case 'S':
        case 'M':
        case 'D':
        case 'N':
        case 'X':
        case '=': pos+= bam_cigar_oplen(cigar[i]);break;
        default: break;
        }
      i++;
      }
  return  scm_from_signed_integer(pos);
	}		
	
static SCM hts_read_cigar(SCM scm_ctx) {
  HtsGuileCtxPtr ptr = cast_to_ctx_ptr(scm_ctx);
  SCM L;
  if(ptr->b->core.n_cigar>0)
    {
    int i=0;
    L = scm_c_make_vector ((int)(ptr->b->core.n_cigar),SCM_UNDEFINED);
    uint32_t *cigar = bam_get_cigar(ptr->b);
    for (i = 0; i < ptr->b->core.n_cigar; ++i) 
        {
        SCM item =  scm_from_pointer(
              (void*)makeCigarElement(
                (int) bam_cigar_oplen(cigar[i]),
                bam_cigar_opchr(cigar[i])
                ),
               free
              );
        scm_vector_set_x (L, scm_from_int (i), item);
        }
    }
  else
    {
    L = scm_c_make_vector (0,SCM_UNDEFINED);
    }
  return L;
  }

static SCM hts_read_cigar_string(SCM scm_ctx) {
    SCM cigar_str;
    HtsGuileCtxPtr ptr = cast_to_ctx_ptr(scm_ctx);
    if (ptr->b->core.n_cigar>0) { // cigar
         int i;
          kstring_t str;
         memset(&str, 0, sizeof(kstring_t));
        uint32_t *cigar = bam_get_cigar(ptr->b);
        for (i = 0; i < ptr->b->core.n_cigar; ++i) {
            kputw(bam_cigar_oplen(cigar[i]), &str);
            kputc(bam_cigar_opchr(cigar[i]), &str);
        	}
        cigar_str =  scm_from_locale_string(str.s);
        free(str.s);
        return cigar_str;
        }
    else
    	{
    	return SCM_UNDEFINED;
    	}
     }

static SCM hts_read_seq_length(SCM scm_ctx) {
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.l_qseq);
	}

static SCM hts_read_seq(SCM scm_ctx) {
	  HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	  if (ptr->b->core.l_qseq) { // seq and qual
	  	kstring_t str;
	  	SCM seq;
	  	int i;
         	memset(&str, 0, sizeof(kstring_t));
		uint8_t *s = bam_get_seq(ptr->b);
		for (i = 0; i < ptr->b->core.l_qseq; ++i) {
			kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], &str);
			}

		seq =  scm_from_locale_string(str.s);		
        	free(str.s);
        	return seq;
	   	}	
	    else
	    	{
	    	return SCM_UNDEFINED;
	    	}
	  }

static SCM hts_read_seq_at(SCM scm_ctx,SCM scm_idx) {
	  HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	  
	  if(!scm_is_integer(scm_idx)) {
	  return SCM_UNDEFINED;
	  }
	  int idx = scm_to_int(scm_idx);
	  if (ptr->b->core.l_qseq && idx < ptr->b->core.l_qseq && idx>=0)  {
		  uint8_t *s = bam_get_seq(ptr->b);
		  char c = "=ACMGRSVTWYHKDBN"[bam_seqi(s, idx)];
		  return   scm_from_char(c);		 	
	   	}	
	    else
	    	{
	    	 fprintf(stderr,"CHAR IS UNDEFINED\n");
	    	return SCM_UNDEFINED;
	    	}
  }

static SCM hts_read_flag(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_signed_integer(ptr->b->core.flag);
	}

static SCM hts_read_is_paired(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool(ptr->b->core.flag & BAM_FPAIRED );
	}

static SCM hts_read_is_proper_pair(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool( (ptr->b->core.flag & BAM_FPAIRED) && (ptr->b->core.flag & BAM_FPROPER_PAIR) );
	}
static SCM hts_read_is_unmapped(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool(ptr->b->core.flag & BAM_FUNMAP );
	}

static SCM hts_mate_is_unmapped(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool( (ptr->b->core.flag & BAM_FPAIRED) && (ptr->b->core.flag & BAM_FMUNMAP) );
	}

static SCM hts_read_reverse_strand(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool(
	  (ptr->b->core.flag & BAM_FREVERSE)
	  );
	}


static SCM hts_mate_reverse_strand(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool(
	  (ptr->b->core.flag & BAM_FPAIRED) && 
	  !(ptr->b->core.flag & BAM_FMUNMAP) && 
	  (ptr->b->core.flag & BAM_FMREVERSE)
	  );
	}

static SCM hts_read_1st_in_pair(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool(
	  (ptr->b->core.flag & BAM_FREAD1)
	  );
	}
static SCM hts_read_2nd_in_pair(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=cast_to_ctx_ptr(scm_ctx);
	return scm_from_bool(
	  (ptr->b->core.flag & BAM_FREAD2)
	  );
	}	
static SCM hts_read_secondary_align(SCM scm_ctx)
	{
	return scm_from_bool(cast_to_ctx_ptr(scm_ctx)->b->core.flag & BAM_FSECONDARY);
	}
static SCM hts_read_qcfail(SCM scm_ctx)
	{
	return scm_from_bool(cast_to_ctx_ptr(scm_ctx)->b->core.flag & BAM_FQCFAIL);
	}
static SCM hts_read_duplicate(SCM scm_ctx)
	{
	return scm_from_bool(cast_to_ctx_ptr(scm_ctx)->b->core.flag & BAM_FDUP);
	}
static SCM hts_read_supplementary(SCM scm_ctx)
	{
	return scm_from_bool(cast_to_ctx_ptr(scm_ctx)->b->core.flag & BAM_FSUPPLEMENTARY);
	}


static void hts_guile_define_module(void *data UNUSED)
	{
	scm_c_define_gsubr ("hts-mate-reverse-strand?", 1, 0, 0, hts_mate_reverse_strand);
	scm_c_define_gsubr ("hts-mate-unmapped?", 1, 0, 0, hts_mate_is_unmapped);
	scm_c_define_gsubr ("hts-read-1st-in-pair?", 1, 0, 0, hts_read_1st_in_pair);
	scm_c_define_gsubr ("hts-read-2nd-in-pair?", 1, 0, 0, hts_read_2nd_in_pair);
  scm_c_define_gsubr ("hts-read-flag", 1, 0, 0, hts_read_flag);
	scm_c_define_gsubr ("hts-read-length", 1, 0, 0, hts_read_seq_length);
	scm_c_define_gsubr ("hts-read-name", 1, 0, 0, hts_read_name);
	scm_c_define_gsubr ("hts-read-pos", 1, 0, 0, hts_read_pos);
	scm_c_define_gsubr ("hts-read-end", 1, 0, 0, hts_read_reference_end);
	scm_c_define_gsubr ("hts-read-unclipped-start", 1, 0, 0, hts_read_unclipped_start);
	scm_c_define_gsubr ("hts-read-unclipped-end", 1, 0, 0, hts_read_unclipped_end);
	scm_c_define_gsubr ("hts-read-reverse-strand?", 1, 0, 0, hts_read_reverse_strand);
	scm_c_define_gsubr ("hts-read-seq", 1, 0, 0, hts_read_seq);
  scm_c_define_gsubr ("hts-read-seq-at", 2, 0, 0, hts_read_seq_at);
  scm_c_define_gsubr ("hts-read-paired?", 1, 0, 0, hts_read_is_paired);
  scm_c_define_gsubr ("hts-read-proper-pair?", 1, 0, 0, hts_read_is_proper_pair);
  scm_c_define_gsubr ("hts-read-unmapped?", 1, 0, 0, hts_read_is_unmapped);

  scm_c_define_gsubr ("hts-read-secondary?", 1, 0, 0, hts_read_secondary_align);
  scm_c_define_gsubr ("hts-read-qcfail?", 1, 0, 0, hts_read_qcfail);
  scm_c_define_gsubr ("hts-read-duplicate?", 1, 0, 0, hts_read_duplicate);
  scm_c_define_gsubr ("hts-read-supplementary?", 1, 0, 0, hts_read_supplementary);
  scm_c_define_gsubr ("hts-read-mapq", 1, 0, 0, hts_read_mapq);
  scm_c_define_gsubr ("hts-read-cigar-string", 1, 0, 0,  hts_read_cigar_string); 
  scm_c_define_gsubr ("hts-read-cigar", 1, 0, 0,  hts_read_cigar); 
  scm_c_define_gsubr ("hts-read-contig", 1, 0, 0,  hts_read_contig); 
  scm_c_define_gsubr ("hts-read-tid", 1, 0, 0,  hts_read_tid); 

	scm_c_export(
	  "hts-mate-reverse-strand?",
	  "hts-mate-unmapped?",
	  "hts-read-1st-in-pair?",
	  "hts-read-2nd-in-pair?",
	  "hts-read-cigar-string",
	  "hts-read-cigar",
	  "hts-read-contig",
	  "hts-read-flag",
		"hts-read-paired?",
		"hts-read-length",
		"hts-read-name",
		"hts-read-pos",
		"hts-read-proper-pair?",
		"hts-read-reverse-strand?",
		"hts-read-seq",
		"hts-read-seq-at",
		"hts-read-tid",
		"hts-read-unmapped?",
		"hts-read-secondary?",
    "hts-read-qcfail?",
    "hts-read-duplicate?",
    "hts-read-supplementary?",
    "hts-read-mapq",
    "hts-read-clipped?",
    "hts-read-unclipped-start",
	  "hts-read-unclipped-end",
	  "hts-read-end",
		NULL);
		
  scm_c_eval_string(__htsguild_scm);
	}


void* hts_guile_init() {
if(scm_is_eq(make_mod,SCM_EOL)) {
		make_mod = scm_c_define_module ("hts", hts_guile_define_module, NULL);
		}
return NULL;
}

static void usage_filter() {
}


static HtsGuileCtxPtr HtsGuileCtxNew()
  {
  HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)malloc(sizeof(HtsGuileCtx));
  memset((void*)ptr,0,sizeof(HtsGuileCtx));
  ptr->b = bam_init1();
  return ptr;
  }

static void* HtsGuileCtxFree(HtsGuileCtxPtr ptr)
  {
  if(ptr==NULL) return NULL;
  if(ptr->b!=NULL) bam_destroy1(ptr->b);
  ptr->b = NULL;
  free(ptr);
  return NULL;
  }

int main_filtersam(int argc,char** argv) {
	char* scriptexpr=NULL;
	char* scriptfile=NULL;
	char* outputfile=NULL;
	int binary=0;
	int compress_level=5;
	samFile *in = 0, *out = 0;

	
	for(;;)
		{
		static struct option long_options[] =
		     {
		      {"expression",  required_argument ,0, 'e'},
			    {"script",   required_argument, 0, 'f'},
			    {"output",   required_argument, 0, 'o'},
			    {"compress", required_argument,0,'c'},
			    {"bam",   no_argument, 0, 'b'},
			    {"help",   no_argument, 0, 'h'},
			    {"version",   no_argument, 0, 'v'},
		      {0, 0, 0, 0}
		     };
		int option_index = 0;
		int c = getopt_long (argc, argv, "hbvf:e:o:c:",
		                    long_options, &option_index);
		if (c == -1) break;
		switch(c)
		  {
		  case 'c' :
		        compress_level = atoi(optarg);
		        if(compress_level<0) compress_level=0;
		        if (compress_level>9) compress_level=9;
		        break;
		  case 'v': printf("version"); return EXIT_SUCCESS;
		  case 'h': usage_filter();return EXIT_SUCCESS;
		  case 'b': binary=1;break;
		  case 'o': outputfile = optarg ; break;
		  case 'e': scriptexpr = optarg ; break;
		  case 'f': scriptfile = optarg ; break;
		  case '?':
		       		fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		       		return EXIT_FAILURE;
	   	default:
	   	        fprintf (stderr, "Bad input.\n");
		       		return EXIT_FAILURE;
		  }
		}
		
		if(scriptexpr==NULL && scriptfile==NULL) {
		  fprintf(stderr,"expression or file must be defined");
		  return EXIT_FAILURE;
		  }
		if(scriptexpr!=NULL && scriptfile!=NULL) {
		  fprintf(stderr,"expression and file both be defined");
		  return EXIT_FAILURE;
		  }

scm_init_guile();
scm_with_guile (hts_guile_init, NULL);
scm_c_use_module("hts");

fprintf(stderr,"compiling %s\n",scriptexpr);
if(scriptexpr!=NULL)
  {
  scm_c_eval_string(scriptexpr);
  }
else
  {
  scm_c_primitive_load(scriptfile);
  }
//scm_c_catch 


SCM filterfun = scm_c_lookup(GUILE_FILTER_NAME);


if(scm_is_false(scm_variable_p(filterfun)))
  {
   fprintf(stderr,"Guile " GUILE_FILTER_NAME " is not a procedure...\n");
	}
SCM filterproc = scm_variable_ref(filterfun);
		htsFormat ga_in = {0,0};
		htsFormat ga_out = {0,0};
		if(optind==argc)
		  {
		  in = hts_open_format("-", "r", &ga_in);
		  }
		else if(optind+1==argc)
	    {
	    in = hts_open_format(argv[optind], "r", &ga_in);
	    }
	  else
	     {
	     fprintf(stderr,"illegal number of args.");
	     return EXIT_FAILURE;
	     }
	 if(in==NULL) {
	  fprintf(stderr,"Cannot open input");
	     return EXIT_FAILURE;
	  }
	  
	  
	  
	  
	  bam_hdr_t* header;
	  if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header");
        return EXIT_FAILURE;
    }
  char output_mode[5];
  int i=0;
  output_mode[i++]='w';
  if(binary==1)
    {
    output_mode[i++]='b';
    output_mode[i++]= compress_level + '0';
    }
  output_mode[i++]='\0';
  
	if ((out = sam_open_format(outputfile?outputfile: "-", output_mode, &ga_out)) == 0)
	  {
	  return EXIT_FAILURE;
	  }
	 
	 if (sam_hdr_write(out, header) != 0) {
                fprintf(stderr, "[main_samview] failed to write the SAM header\n");
              return EXIT_FAILURE;
            }

      int r;
      for(;;)
         {
         HtsGuileCtxPtr ctx = HtsGuileCtxNew();
         ctx->header = header;
         if((r = sam_read1(in, ctx->header, ctx->b)) < 0)
            {
            HtsGuileCtxFree(ctx);
            break;
            }
         SCM sc_ptr = scm_from_pointer((void*)ctx,(scm_t_pointer_finalizer)HtsGuileCtxFree);
         SCM ret=scm_call_1(filterproc,sc_ptr);
         if(scm_is_bool(ret))
            {
      		  if(!scm_to_bool (ret)) continue;
      		  }
        else 
        		{
        		fprintf(stderr,"Guile returned value that is not a boolean\n");
        		continue;
        		}

          r = sam_write1(out, ctx->header, ctx->b);
          if (r < 0)
            {
             fprintf(stderr, "[htsguile] cannot write.\n");
             return EXIT_FAILURE;
             }
            
          if (r < -1) {
             fprintf(stderr, "[htsguile] truncated file.\n");
             return EXIT_FAILURE;
            }
          }
   sam_close(in);
   sam_close(out);
   if ( header ) bam_hdr_destroy(header);
  
return EXIT_SUCCESS;
}

void usage() {

}

int main(int argc,char** argv)
	{
	if(argc<=1) {
	  usage();
	  return EXIT_SUCCESS;	
	  }
	else if(strcmp(argv[1],"filtersam")==0) {
	  return main_filtersam(argc-1,&argv[1]);
	  }
	else 	{
	  fprintf(stderr,"unknown program \"%s\"\n",argv[1]);
	  return EXIT_FAILURE;
    }
	return EXIT_SUCCESS;
	}
