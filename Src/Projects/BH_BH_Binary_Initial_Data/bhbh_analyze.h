#include "bhbh_header.h"
#include "physics_blackhole_lib.h"

/* useful macro */
#define MAX_STR_LEN (400)
#define LINE_STR    "-------------------------------------------------------------------------"


void bhbh_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen);

void bhbh_analyze(Physics_T *const phys,const int iteration);
static void compute_properties(Physics_T *const phys/* bhbh */);



