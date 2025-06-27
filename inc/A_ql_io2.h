int N_ql=15;
static struct A_ql_io2
{
   struct A_proc_info My; /* General IO information */
   int Size;         /* Control: Size of the Shmem */
   int jrho_beg;           /* Input:   */
   int jrho_end;           /* Input:   */
   double* QLarrays;
} *ql_io2;
