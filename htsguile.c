#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <assert.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htsguile.h"

int main_filtersam(int argc,char** argv) {
	char* scriptexpr=NULL;
	char* scriptfile=NULL;
	char* outputfile=NULL;
	
	samFile *in = 0, *out = 0, *un_out=0;

	
	for(;;)
		{
		static struct option long_options[] =
		     {
		      {"expression",  required_argument ,0, 'e'},
			    {"script",   required_argument, 0, 'f'},
			    {"output",   required_argument, 0, 'o'},
			    {"help",   no_argument, 0, 'h'},
			    {"version",   no_argument, 0, 'v'},
		      {0, 0, 0, 0}
		     };
		int option_index = 0;
		int c = getopt_long (argc, argv, "hvf:e:o:",
		                    long_options, &option_index);
		if (c == -1) break;
		switch(c)
		  {
		  case 'v': printf("version"); return EXIT_SUCCESS;
		  case 'h': break;
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
	if ((out = sam_open_format(outputfile?outputfile: "-", "w", &ga_out)) == 0)
	  {
	  return EXIT_FAILURE;
	  }
	 
	 if (sam_hdr_write(out, header) != 0) {
                fprintf(stderr, "[main_samview] failed to write the SAM header\n");
              return EXIT_FAILURE;
            }
    bam1_t *b = bam_init1();
      int r;
      while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
          
          r = sam_write1(out, header, b);
          if (r < 0)
            {
            return EXIT_FAILURE;
             }
            
          if (r < -1) {
            fprintf(stderr, "[main_samview] truncated file.\n");
              return EXIT_FAILURE;
            }
          }
   bam_destroy1(b);
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
