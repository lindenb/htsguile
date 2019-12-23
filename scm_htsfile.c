#include <libguile.h>
#include <htslib/vcf.h>

static void close_file(void* ptr) {
htsFile* f = (htsFile*)ptr;
if(f==NULL) return;
hts_close(f);
}

static SCM scm_hts_open(SCM s_fname,SCM s_mode) {
htsFile* f = NULL;
char* name = NULL;
char* mode = NULL;
SCM_ASSERT (scm_is_string(s_fname), s_fname, 1, "hts-open");
SCM_ASSERT (scm_is_string(s_mode), s_mode, 1, "hts-open");
name = scm_to_locale_string(s_fname);
mode = scm_to_locale_string(s_mode);
f=hts_open(name,mode);
free(name);
free(mode);
return scm_from_pointer(f,close_file);
}

static SCM scm_hts_close(SCM s_file) {
SCM_ASSERT (SCM_POINTER_P(s_file),s_file,1,"hts-close");
if(scm_to_pointer(s_file)==NULL) return SCM_BOOL_F;
close_file(scm_to_pointer(s_file));
return SCM_BOOL_T;
}

void scm_htsfile_init() {
  scm_c_define_gsubr ("hts-open", 2, 0, 0, scm_hts_open);
  scm_c_define_gsubr ("hts-close",1, 0, 0, scm_hts_close);
}
