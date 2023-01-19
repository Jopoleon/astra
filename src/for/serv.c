#include	<unistd.h>
#include	<time.h>
#include	<string.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<unistd.h>
#include	<errno.h>
#include	<sys/times.h>
#include	<sys/types.h>
#include	<sys/ipc.h>
#include	<sys/sem.h>
#include	<sys/shm.h>
#define INT_	int
#define IPCACTIVE 1

void getid (double*,int*,int*);	void	getid_(double*,int*,int*);
void getval (int*,double*,double*);	void	getval_(int*,double*,double*);
void stcopy (char*,char*,int);
void iroundA(char*,int*,int*,int*);
void colovm (int*);

/**************************************************************************/
/* The function returns the total CPU time [sec] spent by the calling process.
      It adds the CPU time spent between two successive calls to the argument.
   The  function "times" returns the number of clock ticks that have elapsed
      since the moment the system was booted. 
   The  "tms_utime"  field contains the CPU time spent executing instructions
      of the calling process.
   The  "tms_stime"  field contains the CPU time spent in the system while 
      executing tasks on behalf of the calling process.                     */
double swatch_(secs)		double	*secs;
{  double	runsec;
   clock_t	  cpu_time, run_time;
   static clock_t time0=0, prev_time;
   static double  secs_per_tick;
   struct   tms	  buf;
   if (time0 == -1) { return -1.; }		/* Overflow range of clock_t*/
   if (time0 ==  0) { 				/* Set time0 at start */
     secs_per_tick = 1./sysconf(_SC_CLK_TCK);
     time0 = times(&buf);
     prev_time = buf.tms_utime+buf.tms_stime;
     return 0.;
   }
   run_time = times(&buf)-time0;		/* Set time difference */
   cpu_time = buf.tms_utime+buf.tms_stime;
   *secs += (cpu_time-prev_time)*secs_per_tick;
   /* printf("Process: %f  %f\n", cpu_time*secs_per_tick, */
   /* (cpu_time-prev_time)*secs_per_tick); */
   prev_time = cpu_time;	
   runsec = run_time*secs_per_tick;   /* printf("Process: %f\n", runsec); */
   return runsec;
}
double swatch (secs)		double	*secs;
{  double	runsec;
   runsec = swatch_(secs);
   return runsec;
}

/***************** Get ID of calling function and its parent  *************/
void getid (var, mediator, id)	double *var;	int *mediator, *id;
{    getid_(var, mediator, id);	}
void getid_(var, mediator, id)	double *var;	int *mediator, *id;
{	/* No more than 16 functions for no more than 2200 calling each */
	static	double	*varid[16][2200];
	static	int	*callid[16];
	int		i, j;
	if   ( *mediator == 0 )
	     { for ( i=0; i<=15; ) { if ( callid[i] == 0) break; i++; }
		if ( i==16 )	goto Message1;
		callid[i] = mediator;	/* identify & store calling function */
		*mediator=i+1;
	     }				/* loop over calling functions */
	if   ( *mediator < 0 || *mediator >= 17)	goto Message2;

	i = *mediator-1;
	for (j=0; j<=2109;) { if (varid[i][j] == var) { *id=++j; return; }
	if ( varid[i][j] == 0  ) break; j++;
			  }	if ( j==2200 )	goto Message3;
	varid[i][j] = var; *id = ++j;
	return;

Message1:  printf(">>> GetID:  > 16 calling functions\n");	  return;
Message2:  printf(">>> GetID:  Illegal 2nd parameter\n"); *id=-2; return;
Message3:  /* Too many varibles > 2200 */		  *id=-1; return;
}
/****************** Get ID of calling function and its parent  *************/
void getval (shift, arr, value)	double *value, *arr;	int *shift;
{    getval_(shift, arr, value);	}
void getval_(shift, arr, value)	double *value, *arr;	int *shift;
{	*value = *(arr+*shift);
	/* printf("Shift = %d,    -> = %lf\n",*shift,*value); */
	return;
}
/**************************************************************************/
void stcopy (s1,s2,n) int n; char *s1, *s2;
/* Copy n characters from string s2 to string s1 
   equivalent to	{ strncpy(s1, s2, n);	s1[n] = '\0'; }  */
{	int	i;
	for(i=0; i <= n-1 ; ++i)	*(s1+i) = *(s2+i);
	*(s1+i) = '\0';
}
/**************************************************************************/
int numstrA(x, str, l)	char *str; double x; int l;
/* Inpit:
	x   - number
	l   - string length
  Output:
	str - string of the length "l"
*/
{	static int      n = 12;
	int             i, ne, is, ll, k, i5;
	char            ch[30];

	for (i = 0; i < l; i++) str[i] = ' ';
	str[l] = '\0';
	if (x == 0.) { str[l - 1] = '0'; return (0); }
	sprintf(ch, "%+.*e", n, x);
	i = sscanf((ch + n + 4), "%d", &ne);
	if (ch[0] == '-') is = 1;
	else		  is = 0;
	ll = l - is;
	if (ll <= 0) { str[l - 1] = '*'; return (1); }
	k = n + 2;
	while (ch[k] == '0') k--;
	ch[0] = ch[1];
	for (i = 1; i < k; ch[i] = ch[i + 2], i++); k--;
	ch[k] = '\0'; ne += (1 - k); i5 = 0;
	if (k > ll) { ne = ne + k - ll; k = ll; iroundA(ch, &k, &ne, &i5); }
	while (k) {
	    if (ne >= 0) {
		if (ne + k <= ll) {
		    for (i = l - 1; ne > 0; ne--, str[i--] = '0');
		    for (k--; k >= 0; str[i--] = ch[k--]);
		    if (is) str[i] = '-'; return (0); } 
		else {	i = 2; if (ne > 9) i++; if (ne > 99) i++;
		    if (k + i > ll) { k--;
			if (k==0) { str[l-1]='*'; return (1); }
			ne++; iroundA(ch, &k, &ne, &i5); } 
		    else { i = l - i; sprintf(str + i, "e%d", ne);
			for (--k; k >= 0; str[--i]=ch[k--]);
			if (is) str[--i] = '-';
			return (0); } } } 
	    else { if (ll > -ne) {
		    if (k - ll) { i = k + ne;
			if (i >= 0) { i = l; while (ne) 
			    { k--; i--; str[i] = ch[k]; ne++; }
			    i--; str[i] = '.';
			    while (k) { k--; i--; str[i] = ch[k]; } } 
			else { i = l; while (k) 
			    { k--; i--; str[i]=ch[k]; ne++; }
			    while (ne) { i--; str[i]='0'; ne++; }
			    i--; str[i] = '.'; }
			if (is) { i--; str[i] = '-'; }
			return (0); } 
		    else { k--; ne++; iroundA(ch, &k, &ne, &i5); } }
		else {	i = 3; if (ne < -9) i = 4; if (ne < -99) i = 5;
		    if (ll >= k + i) { i = l - i;
			sprintf(str + i, "e%d", ne);
			while (k) { i--; k--; str[i] = ch[k]; }
			if (is) { i--; str[i] = '-'; }
			return (0); }
		    k--; if (k == 0) { str[l-1] = '*'; return (1); }
		    ne++; iroundA(ch, &k, &ne, &i5);
	}   }	}
	return (0);
}
/**************************************************************************/
int num2str_(x1, str, l1)	char *str; double *x1; int *l1;
/* the same as numstrA but with float argument and callable from FORTRAN */
{	static int      n = 12;
	int             l, i, ne, is, ll, k, i5;
	char            ch[30];
	double		x;

	l = *l1;	x = *x1;
	for (i = 0; i < l; i++) str[i] = ' ';
	str[l] = '\0';
	if (x == 0.) { str[l - 1] = '0'; return (0); }
	sprintf(ch, "%+.*e", n, x);
	i = sscanf((ch + n + 4), "%d", &ne);
	if (ch[0] == '-') is = 1;
	else		  is = 0;
	ll = l - is;
	if (ll <= 0) { str[l - 1] = '*'; return (1); }
	k = n + 2;
	while (ch[k] == '0') k--;
	ch[0] = ch[1];
	for (i = 1; i < k; ch[i] = ch[i + 2], i++); k--;
	ch[k] = '\0'; ne += (1 - k); i5 = 0;
	if (k > ll) { ne = ne + k - ll; k = ll; iroundA(ch, &k, &ne, &i5); }
	while (k) {
	    if (ne >= 0) {
		if (ne + k <= ll) {
		    for (i = l - 1; ne > 0; ne--, str[i--] = '0');
		    for (k--; k >= 0; str[i--] = ch[k--]);
		    if (is) str[i] = '-'; return (0); } 
		else {	i = 2; if (ne > 9) i++; if (ne > 99) i++;
		    if (k + i > ll) { k--;
			if (k==0) { str[l-1]='*'; return (1); }
			ne++; iroundA(ch, &k, &ne, &i5); } 
		    else { i = l - i; sprintf(str + i, "e%d", ne);
			for (--k; k >= 0; str[--i]=ch[k--]);
			if (is) str[--i] = '-';
			return (0); } } } 
	    else { if (ll > -ne) {
		    if (k - ll) { i = k + ne;
			if (i >= 0) { i = l; while (ne) 
			    { k--; i--; str[i] = ch[k]; ne++; }
			    i--; str[i] = '.';
			    while (k) { k--; i--; str[i] = ch[k]; } } 
			else { i = l; while (k) 
			    { k--; i--; str[i]=ch[k]; ne++; }
			    while (ne) { i--; str[i]='0'; ne++; }
			    i--; str[i] = '.'; }
			if (is) { i--; str[i] = '-'; }
			return (0); } 
		    else { k--; ne++; iroundA(ch, &k, &ne, &i5); } }
		else {	i = 3; if (ne < -9) i = 4; if (ne < -99) i = 5;
		    if (ll >= k + i) { i = l - i;
			sprintf(str + i, "e%d", ne);
			while (k) { i--; k--; str[i] = ch[k]; }
			if (is) { i--; str[i] = '-'; }
			return (0); }
		    k--; if (k == 0) { str[l-1] = '*'; return (1); }
		    ne++; iroundA(ch, &k, &ne, &i5);
	}   }	}
	return (0);
}
/**************************************************************************/
void iroundA(ch, k, ne, i5)	char *ch;	int *k, *ne, *i5;
{	int             i;
	i = *k - 1;
	if (ch[*k] >= 53 + *i5) { ch[i]++; *i5 = 0; if (ch[i] == 53) *i5 = 1; }
	while (ch[i] == 58) {
		if (i > 0) { ch[i] = 48; i--; ch[i]++; *i5 = 0;
			if (ch[i] == 53) *i5 = 1;
			(*k)--; (*ne)++; } 
		else { ch[i] = 49; (*ne)++; } }
	while (ch[(*k) - 1] == 48) { if ((*k) > 1) { (*k)--; (*ne)++; } }
}


void fenvex(void);			void fenvex_(void);
void fenvdx(void);			void fenvdx_(void);
#ifdef AFENV
#include <signal.h>
#include <fenv.h>
void float_error(int);
void sigint_handler(int sig); /* prototype */
void fenvex_(void)
{
  feenableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO );		  /*
  feenableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); */
  /* @ Linux FE_INVALID=1,  FE_OVERFLOW=8,  FE_DIVBYZERO=4,  FE_UNDERFLOW=16
  printf("%d %d %d %d\n",FE_INVALID,FE_OVERFLOW,FE_DIVBYZERO,FE_UNDERFLOW);
  */
  signal(SIGFPE, float_error);
  signal(SIGINT, sigint_handler);
  return;
}
void fenvex(void)
{ void fenvex_(); }

void fenvdx_(void)
{
  fedisableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO );		  /*
  fedisableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); */
  signal(SIGFPE, float_error);
  return;
}
void fenvdx(void)
{ void fenvdx_(); }

void sigint_handler(int sig)			/* this is the handler */
{
  const int j = 47;  
  /*  printf("Not this time!\n"); */
  /* system("if ( -e ${HOME}/bin/cln ) ${HOME}/bin/cln"); */
  /* if (system("/../ESC/cln") != -1)  ifkey_(&j); */
  ifkey_(&j);
}

void float_error(int i)
{
  int	j;
  fprintf(stderr,"\n >>> Floating point exception in\n");
  j = 257;	ifkey_(&j);
  exit(EXIT_FAILURE);
}
#else
void fenvex_(void)
{  return;	  }
void fenvex(void)
{ void fenvex_(); }
void fenvdx_(void)
{  return;	  }
void fenvdx(void)
{ void fenvdx_(); }
#endif

