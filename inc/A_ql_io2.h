static struct A_ql_io2
{
   struct A_proc_info My; /* General IO information */
   int Size;         /* Control: Size of the Shmem */
   int is;           /* Input:   */
   int ie;           /* Input:   */
   double *chi;  /* Output: chi_i  */
   double *che;  /* Output: chi_e  */
   double *dif;  /* Output: diff_i  */
   double *vin;  /* Output: diff_i  */
   double *dph;  /* Output: diff_tor  */
   double *dpl;  /* Output: diff_pll  */
   double *dpr;  /* Output: diff_perp  */
   double *xtb;  /* Output: turb heat exchange  */
   double *egm;  /* Output: egamma  */
   double *gam;  /* Output: gamma_p  */
   double *gm1;  /* Output: lead mode rate  */
   double *gm2;  /* Output: 2nd mode rate  */
   double *om1;  /* Output: lead mode freq  */
   double *om2;  /* Output: 2nd mode freq  */
   double *fr1;  /* Output:   */
} *IOQL;
