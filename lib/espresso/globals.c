/*
 * Revision Control Information
 *
 * $Source$
 * $Author$
 * $Revision$
 * $Date$
 *
 */
#include "espresso.h"




/*
 *    Global Variable Declarations
 */

unsigned int debug;              /* debug parameter */
BOOL verbose_debug;              /* -v:  whether to print a lot */
char *total_name[TIME_COUNT];    /* basic function names */
long total_time[TIME_COUNT];     /* time spent in basic fcts */
int total_calls[TIME_COUNT];     /* # calls to each fct */

BOOL echo_comments;         /* turned off by -eat option */
BOOL echo_unknown_commands;     /* always true ?? */
BOOL force_irredundant;          /* -nirr command line option */
BOOL skip_make_sparse;
BOOL kiss;                       /* -kiss command line option */
BOOL pos;                        /* -pos command line option */
BOOL print_solution;             /* -x command line option */
BOOL recompute_onset;            /* -onset command line option */
BOOL remove_essential;           /* -ness command line option */
BOOL single_expand;              /* -fast command line option */
BOOL summary;                    /* -s command line option */
BOOL trace;                      /* -t command line option */
BOOL unwrap_onset;               /* -nunwrap command line option */
BOOL use_random_order;         /* -random command line option */
BOOL use_super_gasp;         /* -strong command line option */
char *filename;             /* filename PLA was read from */

struct pla_types_struct pla_types[] = {
    {"-f", F_type},
    {"-r", R_type},
    {"-d", D_type},
    {"-fd", FD_type},
    {"-fr", FR_type},
    {"-dr", DR_type},
    {"-fdr", FDR_type},
    {"-fc", F_type | CONSTRAINTS_type},
    {"-rc", R_type | CONSTRAINTS_type},
    {"-dc", D_type | CONSTRAINTS_type},
    {"-fdc", FD_type | CONSTRAINTS_type},
    {"-frc", FR_type | CONSTRAINTS_type},
    {"-drc", DR_type | CONSTRAINTS_type},
    {"-fdrc", FDR_type | CONSTRAINTS_type},
    {"-pleasure", PLEASURE_type},
    {"-eqn", EQNTOTT_type},
    {"-eqntott", EQNTOTT_type},
    {"-kiss", KISS_type},
    {"-cons", CONSTRAINTS_type},
    {"-scons", SYMBOLIC_CONSTRAINTS_type},
    {0, 0}
};


struct cube_struct cube, temp_cube_save;
struct cdata_struct cdata, temp_cdata_save;

int bit_count[256] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};


