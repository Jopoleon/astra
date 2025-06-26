static struct A_ql_io2
{
   struct A_proc_info My; /* General IO information */
   int Size;         /* Control: Size of the Shmem */
   int is;           /* Input:   */
   int ie;           /* Input:   */
   double* QLarrays;
} *ql_io2;
