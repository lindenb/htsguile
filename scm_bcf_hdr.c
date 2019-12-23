#include <libguile.h>
#include <htslib/vcf.h>

static void dispose(void* ptr) {
if(ptr==NULL) return;
bcf_hdr_destroy((bcf_hdr_t*)ptr);
}

static SCM scm_bcf_hdr_samples(SCM s_h) {
bcf_hdr_t* h = NULL;
int i=0;
SCM list;
SCM_ASSERT (SCM_POINTER_P(s_h),s_h,2,"bcf-hdr-samples");
h=(bcf_hdr_t*)scm_to_pointer(s_h);
list = scm_c_make_vector((size_t)bcf_hdr_nsamples(h),SCM_UNDEFINED);
for(i=0;i< bcf_hdr_nsamples(h); i++) {
	 scm_c_vector_set_x(list,(size_t)i,scm_from_locale_string(h->samples[i]));
	}
return list;
}

static SCM scm_bcf_hdr_sample(SCM s_h,SCM s_idx) {
bcf_hdr_t* h = NULL;
int i=0;
SCM list;
SCM_ASSERT (SCM_POINTER_P(s_h),s_h,1,"bcf-hdr-sample");
SCM_ASSERT (scm_is_integer(s_idx),s_idx,2,"bcf-hdr-sample");
h=(bcf_hdr_t*)scm_to_pointer(s_h);
i = scm_to_int(s_idx);
if(i<0 || i>=bcf_hdr_nsamples(h)) return SCM_UNDEFINED;
return scm_from_locale_string(h->samples[i]);
}

static SCM scm_bcf_hdr_nsamples(SCM s_h) {
SCM_ASSERT (SCM_POINTER_P(s_h),s_h,1,"bcf-hdr-nsamples");
return scm_from_int(bcf_hdr_nsamples((bcf_hdr_t*)scm_to_pointer(s_h)));
}


static SCM scm_bcf_hdr_write(SCM s_fout,SCM s_h) {
bcf_hdr_t* h = NULL;
htsFile* fp = NULL;
SCM_ASSERT (SCM_POINTER_P(s_fout),s_fout,1,"hts-hdr-write");
SCM_ASSERT (SCM_POINTER_P(s_fout),s_h,2,"hts-hdr-write");

fp=(htsFile*)scm_to_pointer(s_fout);
h=(bcf_hdr_t*)scm_to_pointer(s_h);

return  bcf_hdr_write(fp,h) == 0 ? SCM_BOOL_T:SCM_BOOL_F;
}

static SCM scm_bcf_read(SCM s_fin,SCM s_h) {
htsFile* fp = NULL;
bcf_hdr_t* h = NULL;
bcf1_t* rec = NULL;
SCM_ASSERT (SCM_POINTER_P(s_fin),s_fin,1,"bcf-read");
SCM_ASSERT (SCM_POINTER_P(s_h),s_h,2,"bcf-read");
fp=(htsFile*)scm_to_pointer(s_fin);
h=(bcf_hdr_t*)scm_to_pointer(s_h);
rec = bcf_init();
if(bcf_read(fp,h,rec)!=0) return SCM_UNDEFINED;
return scm_from_pointer(rec,(scm_t_pointer_finalizer)bcf_destroy);
}

static SCM scm_bcf_hdr_read(SCM s_fin) {
bcf_hdr_t* h = NULL;
SCM_ASSERT (SCM_POINTER_P(s_fin),s_fin,1,"bcf-hdr-read");
htsFile* fp=(htsFile*)scm_to_pointer(s_fin);
h = bcf_hdr_read(fp);
return scm_from_pointer(h,dispose);
}

static SCM scm_bcf_header_init(SCM s_mode) {
bcf_hdr_t* h = NULL;
char* mode = NULL;
SCM_ASSERT (scm_is_string(s_mode), s_mode, 1, "bcf-hdr-init");
mode = scm_to_locale_string(s_mode);
h = bcf_hdr_init(mode);
free(mode);
return scm_from_pointer(h,dispose);
}

static SCM scm_bcf_hdr_destroy(SCM s_h) {
dispose(scm_to_pointer(s_h));
return SCM_BOOL_T;
}

void scm_bcf_hdr_init() {
  scm_c_define_gsubr ("bcf-read", 2, 0, 0,scm_bcf_read);
  scm_c_define_gsubr ("bcf-hdr-nsamples", 1, 0, 0,scm_bcf_hdr_samples);
  scm_c_define_gsubr ("bcf-hdr-samples", 1, 0, 0,scm_bcf_hdr_samples);
  scm_c_define_gsubr ("bcf-hdr-sample", 2, 0, 0,scm_bcf_hdr_sample);
  scm_c_define_gsubr ("bcf-hdr-write", 2, 0, 0,scm_bcf_hdr_write);
  scm_c_define_gsubr ("bcf-hdr-read", 1, 0, 0,scm_bcf_hdr_read);
  scm_c_define_gsubr ("bcf-hdr-init", 1, 0, 0,scm_bcf_header_init);
  scm_c_define_gsubr ("bcf-hdr-destroy",1, 0, 0,scm_bcf_hdr_destroy);
}
