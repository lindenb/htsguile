#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include "bcfguile.h"

int main(int argc,char** argv) {
int c;
char* filenameout = NULL;
char* guilexpr = NULL;
struct VariantContext ctx;

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
	if ( bcf_write(vcfout,ctx.header,ctx.rec)!=0 ) {
		fprintf(stderr,"I/O error cannot write record.\n");
		return EXIT_FAILURE;
		}
	}
bcf_hdr_destroy(ctx.header);
bcf_destroy(ctx.rec);
hts_close(vcfin);
hts_close(vcfout);
return EXIT_SUCCESS;
}
