#include <libguile.h>
extern void scm_htsfile_init();
extern void scm_bcf_hdr_init();
extern void scm_bcf_rec_init();
static void inner_main (void *closure, int argc, char **argv)
{
 scm_bcf_rec_init();
 scm_htsfile_init();
 scm_bcf_hdr_init();
  scm_shell (argc, argv);
}

int
main (int argc, char **argv)
{
  scm_boot_guile (argc, argv, inner_main, 0);
  return 0; /* never reached */
}
