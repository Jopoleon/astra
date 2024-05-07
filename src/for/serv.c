#include <unistd.h>
#include <stdio.h>
#include <sys/times.h>

void getid_(double*, int*, int*);
void stcopy(char*, char*, int);
void iroundA(char*, int*, int*, int*);

/**************************************************************************/
/* The function returns the total CPU time [sec] spent by the calling process.
      It adds the CPU time spent between two successive calls to the argument.
   The  function "times" returns the number of clock ticks that have elapsed
      since the moment the system was booted. 
   The  "tms_utime"  field contains the CPU time spent executing instructions
      of the calling process.
   The  "tms_stime"  field contains the CPU time spent in the system while 
      executing tasks on behalf of the calling process.                     */
double swatch_(double *secs){
    double runsec;
    clock_t cpu_time, run_time;
    static clock_t time0=0, prev_time;
    static double secs_per_tick;
    struct tms buf;
    if (time0 == -1){
        return -1.;
    }  /* Overflow range of clock_t*/
    if (time0 ==  0){     /* Set time0 at start */
        secs_per_tick = 1./sysconf(_SC_CLK_TCK);
        time0 = times(&buf);
        prev_time = buf.tms_utime+buf.tms_stime;
        return 0.;
    }
    run_time = times(&buf) - time0;  /* Set time difference */
    cpu_time = buf.tms_utime + buf.tms_stime;
    *secs += (cpu_time - prev_time)*secs_per_tick;
    prev_time = cpu_time; 
    runsec = run_time*secs_per_tick;
    return runsec;
}

/**************************************************************************/
void getid_(double *var, int *mediator, int *id){
/* Get ID of calling function and its parent
   No more than 16 functions for no more than 2200 calling each */
    static double *varid[16][2200];
    static int *callid[16];
    int i, j;

    if (*mediator == 0){
        for (i=0; i<16; i++){
            if (callid[i] == 0) break;
        }
        if (i == 16){
            printf(">>> GetID:  > 16 calling functions\n");
            return;
        }
        callid[i] = mediator; /* identify & store calling function */
        *mediator = i + 1;
    }    /* loop over calling functions */
    if (*mediator < 0 || *mediator >= 17){
        printf(">>> GetID:  Illegal 2nd parameter\n");
        *id = -2;
        return;
    }
    i = *mediator - 1;
    for (j=0; j<2200; j++){
        if (varid[i][j] == var){
            *id = ++j;
            return;
        }
        if (varid[i][j] == 0) break;
    }
    if (j == 2200){
        printf("Too many varibles > 2200\n");
        *id = -1;
    }
    else{
        varid[i][j] = var;
        *id = ++j;
    }
    return;
}

/**************************************************************************/
void stcopy(char *s1, char *s2, int n){
/* Copy n characters from string s2 to string s1 
   equivalent to { strncpy(s1, s2, n); s1[n] = '\0'; }  */
    int i;
    for (i=0; i<n; ++i){
        *(s1+i) = *(s2+i);
    }
    *(s1+i) = '\0';
}

/**************************************************************************/
int num2str(double x, char *str, int l){
/* 
Input:
   x : number
   l : string length
Output:
   str : string of length "l"
*/

    static int n=12;
    int i, ne, is, ll, k, i5;
    char ch[30];

    for (i = 0; i < l; i++){
        str[i] = ' ';
    }
    str[l] = '\0';
    if (x == 0.){
        str[l - 1] = '0';
        return 0;
    }
    sprintf(ch, "%+.*e", n, x);
    i = sscanf((ch + n + 4), "%d", &ne);
    if (ch[0] == '-'){
        is = 1;
    }
    else{
        is = 0;
    }
    ll = l - is;
    if (ll < 1){
        str[l-1] = '*';
        return 1;
    }
    k = n + 2;
    while (ch[k] == '0'){
        k--;
    }
   ch[0] = ch[1];
   for (i=1; i < k; ch[i]=ch[i+2], i++){
       k--;
    }
    ch[k] = '\0'; ne += (1 - k);
    i5 = 0;
    if (k > ll){
        ne = ne + k - ll;
        k = ll;
        iroundA(ch, &k, &ne, &i5);
    }

    while (k) {
        if (ne >= 0) {
            if (ne + k <= ll) {
                for (i=l-1; ne > 0; ne--, str[i--] = '0');
                for (k--; k >= 0; str[i--] = ch[k--]);
                if (is) str[i] = '-';
                return 0;
            }
            else{
                i = 2;
                if (ne >  9) i++;
                if (ne > 99) i++;
                if (k + i > ll){
                    k--;
                    if (k == 0){
                        str[l-1] = '*';
                        return 1;
                    }
                    ne++;
                    iroundA(ch, &k, &ne, &i5);
                }
                else{
                    i = l - i;
                    sprintf(str + i, "e%d", ne);
                    for (--k; k >= 0; str[--i]=ch[k--]);
                    if (is) str[--i] = '-';
                    return 0;
                }
            }
        } 
        else{
            if (ll > -ne){
                if (k - ll){
                    i = k + ne;
                    if (i >= 0){
                        i = l;
                        while (ne){
                            k--;
                            i--;
                            str[i] = ch[k];
                            ne++;
                        }
                        i--;
                        str[i] = '.';
                        while (k){
                            k--;
                            i--;
                            str[i] = ch[k];
                        }
                    }
                    else{
                        i = l;
                        while (k){
                            k--;
                            i--;
                            str[i] = ch[k];
                            ne++;
                        }
                        while (ne){
                            i--;
                            str[i] = '0';
                            ne++;
                        }
                        i--;
                        str[i] = '.';
                    }
                    if (is){
                        i--;
                        str[i] = '-';
                    }
                    return 0;
                }
                else{
                    k--;
                    ne++;
                    iroundA(ch, &k, &ne, &i5);
                }
            }
            else{
                i = 3;
                if (ne <  -9) i = 4;
                if (ne < -99) i = 5;
                if (ll >= k + i){
                    i = l - i;
                    sprintf(str + i, "e%d", ne);
                    while (k){
                        i--;
                        k--;
                        str[i] = ch[k];
                    }
                    if (is){
                        i--;
                        str[i] = '-';
                    }
                    return 0;
                }
                k--;
                if (k == 0){
                    str[l-1] = '*';
                    return 1;
                }
                ne++;
                iroundA(ch, &k, &ne, &i5);
            }
        }
    }
    return 0;
}

/**************************************************************************/
int num2str_(double *x1, char *str, int *l1){
    return num2str(*x1, str, *l1);
}

/**************************************************************************/
void iroundA(char *ch, int *k, int *ne, int *i5){
    int i;
    i = *k - 1;
    if (ch[*k] >= 53 + *i5){
        ch[i]++;
        *i5 = 0;
        if (ch[i] == 53) *i5 = 1;
    }
    while (ch[i] == 58){
        if (i > 0){
            ch[i] = 48;
            i--;
            ch[i]++;
            *i5 = 0;
            if (ch[i] == 53) *i5 = 1;
            (*k)--;
            (*ne)++;
        }
        else{
            ch[i] = 49;
            (*ne)++;
        }
    }
    while (ch[(*k) - 1] == 48){
        if ((*k) > 1){
            (*k)--;
            (*ne)++;
        }
    }
}


void fenvex_(void);

#ifdef AFENV

#include <signal.h>
#include <fenv.h>
void float_error(int);
void sigint_handler(int sig); /* prototype */

void fenvex_(void){
    feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
/*
    feenableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);
    @ Linux FE_INVALID=1,  FE_OVERFLOW=8,  FE_DIVBYZERO=4,  FE_UNDERFLOW=16
    printf("%d %d %d %d\n",FE_INVALID,FE_OVERFLOW,FE_DIVBYZERO,FE_UNDERFLOW);
*/
    signal(SIGFPE, float_error);
    signal(SIGINT, sigint_handler);
    return;
}

void sigint_handler(int sig){   /* this is the handler */
    const int j = 47;  
    ifkey_(&j);
}

void float_error(int i){
    int j;
    fprintf(stderr,"\n >>> Floating point exception in\n");
    j = 257;
    ifkey_(&j);
    exit(EXIT_FAILURE);
}

#else

void fenvex_(void){
    return;
}

#endif
