static struct A_vars
{   /* Common block A_VARIABLES */
  double ab;   		/* AB   */
  double abc;  		/* ABC  */
  double aim1;   	/* AIM1 */
  double aim2;   	/* AIM2 */
  double aim3;   	/* AIM3 */
  double amj;    	/* AMJ,   AWALL*/
  double btor;   	/* BTOR */
  double elong;  	/* ELONG, ELONM*/
  double encl;   	/* ENCL */
  double enwm;   	/* ENWM, FECR,FFW,FICR,FLH,GN2E,GN2I */
  double ipl;    	/* IPL,  LEXT*/
  double nncl;   	/* NNCL */
  double nnwm;   	/* NNWM, QECR,QFW,QICR,QLH,QNBI*/
  double rtor;   	/* RTOR */
  double shift;  	/* SHIFT*/
  double trian;  	/* TRIAN, TRICH, UEXT*/
  double updwn;  	/* UPDWN, WNE,WTE,WTI*/
  double zmj;	  	/* ZMJ  */
	/* Common block A_GRIDS */
  double time;  	/* TIME */
  double hro;   	/* HRO  */
  double roc;   	/* ROC  */
  double tau;   	/* TAU  */	   
  int na1; 	  	/* NA1 */ 
  int nab;	  	/* NAB */
  int nb1;	  	/* NB1 */
  int nrd;	  	/* NRD */
  int na1n;
  int na1e;
  int na1i; 
  int	 c_size[12];  /* sizeof: int, double, key_t, pid_t, etc */
  int    CheckWord;     /* Control number */
} *AVARS;
