#include <unistd.h>
#include <stdio.h>
#include <sys/times.h>

void getid_(double*, int*, int*);
void fenvex_(void);

/**********************************************************************/
double swatch_(double *secs){
/* The function returns the total CPU time [sec] spent by the calling process.
      It adds the CPU time spent between two successive calls to the argument.
   The  function "times" returns the number of clock ticks that have elapsed
      since the moment the system was booted. 
   The  "tms_utime"  field contains the CPU time spent executing instructions
      of the calling process.
   The  "tms_stime"  field contains the CPU time spent in the system while 
      executing tasks on behalf of the calling process.               */
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

/**********************************************************************/
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


#ifdef AFENV

#include <signal.h>
#include <fenv.h>
void float_error(int);
void sigint_handler(int); /* prototype */

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
