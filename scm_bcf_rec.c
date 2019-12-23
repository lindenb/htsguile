#include <libguile.h>
#include <htslib/vcf.h>

static void dispose(void* ptr) {
if(ptr==NULL) return;
bcf_destroy((bcf1_t*)ptr);
}

static SCM scm_bcf_record_init() {
bcf1_t* b = bcf_init();
return scm_from_pointer(b,dispose);
}

static SCM scm_bcf_rec_destroy(SCM s_b) {
dispose(scm_to_pointer(s_b));
return SCM_BOOL_T;
}

void scm_bcf_rec_init() {
  scm_c_define_gsubr ("bcf-init", 0, 0, 0,scm_bcf_record_init);
  scm_c_define_gsubr ("bcf-destroy", 0, 0, 0,scm_bcf_rec_destroy);
}
