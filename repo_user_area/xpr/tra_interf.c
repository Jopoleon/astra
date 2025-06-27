/*---------------------------------------------------------------------*/
/* The file includes functions (double)`swatch' and
   the starting part of the function (int)`sbp2shm' that
   acquires information about AWD, process wd, id, name, key,
   ord. number, allocated shared memory ID and semaphore ID.        */
/*------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <errno.h>
#include "A_vars.h"
#include "A_proc.h"

#ifndef INT8
#define INT_ int
#else
#define INT_ long
#endif

union semun {
    int val;                    /* value for SETVAL */
    struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
    unsigned short int *array;  /* array for GETALL, SETALL */
    struct seminfo *__buf;      /* buffer for IPC_INFO */
};

void a_stop_();
double swatch (double*);
extern INT_ A_NA1;
int   SemID ,  ShMid0,  ShMid1;
void *ShmAd0, *ShmAd1, *ShmAdr;
void qlk_interf_();
void neo_interf_();
void tglf_interf_();
void write_aipc();

#define NC1 A_NA1

/*
  Returns CPU time [sec] between two successive calls to the argument.
  The  function "times" returns the number of clock ticks that have elapsed
     since the moment the system was booted.
  The  "tms_utime"  field contains the CPU time spent executing instructions
     of the calling process.
  The  "tms_stime"  field contains the CPU time spent in the system while
     executing tasks on behalf of the calling process.
*/

double swatch(double *secs){
    double runsec;
    clock_t cpu_time, run_time;
    static clock_t time0=0, prev_time;
    static double secs_per_tick;
    struct tms buf;
    if (time0 == -1){
        return -1.;
    }  /* Overflow range of clock_t*/
    if (time0 == 0){     /* Set time0 at start */
        secs_per_tick = 1./sysconf(_SC_CLK_TCK);
        time0 = times(&buf);
        prev_time = buf.tms_utime + buf.tms_stime;
        return 0.;
    }
    run_time = times(&buf) - time0;  /* Set time difference */
    cpu_time = buf.tms_utime + buf.tms_stime;
    *secs += (cpu_time - prev_time)*secs_per_tick;
    prev_time = cpu_time;
    runsec = run_time*secs_per_tick;
    return runsec;
}

/*---------------------------------------------------------------------*/
/* Note description of arg2, arg3, arg4 should be compatible with
    INT_ or int
*/

int sbp2shm_(char* arg0, char* arg1, int* arg2, int* arg3, int* arg4)
{
  int i, qlSize;
    static union semun Mysemun;
// sembuf members: {sem_num,sem_op,sem_flag};
    static struct sembuf buf0 = {0, 1, IPC_NOWAIT};
    static struct sembuf bufN = {1,-1, ~SEM_UNDO&~IPC_NOWAIT};
    static char whoami[32];
    static char AWD[96];
    static struct A_proc_info Mama, My;

    swatch(&(My.CPUse));
/* Analyze the calling command string. Get own PID and name. */
    My.Pid = getpid();
    getcwd(My.Path, (size_t)64);
    strcat(My.Path, "/");

    strcat(My.Path, arg0+2);
    if (strrchr(arg0, '/') != NULL){
        sscanf(strrchr(arg0, '/')+1, "%s", whoami); }
    else{
        sscanf(arg0, "%s", whoami);
    }
    sscanf(arg1, "%s", Mama.Path);
    Mama.Pid = (pid_t)*arg2;
    Mama.Key = (key_t)*arg3;
    My.OrdNr = *arg4;
    printf("Fortran main: %s %d %d\n", My.Path, Mama.Pid, Mama.Key);

/* Associate My semaphore with the ordinal process number */
    bufN.sem_num = My.OrdNr;
/* Get semaphore and shmem IDs. */
    SemID  = semget(Mama.Key, 0, 0660);
    ShMid0 = shmget((key_t)(Mama.Key+0), 0, 0660);
    ShMid1 = shmget((key_t)(Mama.Key+1), 0, 0660);
    ShmAd0 = shmat(ShMid0, NULL, 0);
    ShmAd1 = shmat(ShMid1, NULL, 0);
    AVARS = (struct A_vars *)ShmAd0;
    A_NA1 = AVARS->na1;
#include "A_arrs.h"
#include "A_ql_IO.h"
    strcpy(AWD, Mama.Path);
    if (strstr(AWD, "bin/") == NULL){
        printf("SBP launch string error\n");
        a_stop_();
    }
    else{
        *strstr(AWD, "bin/") = '\0';
    }
    My.Key = ftok( My.Path, (int)My.Pid);
#include "A_ql_io2.h"
    qlSize = sizeof(struct A_ql_io2) - sizeof(double) + N_ql*NC1*sizeof(double);

/*------------------------------------
  Create My shared memory segment
*/
    My.ShMid = shmget(My.Key, qlSize, 0660|IPC_CREAT);
/* Attach My shared memory to the process
   My shared memory segment starts at ShmAdr */
    ShmAdr = shmat(My.ShMid, NULL, 0);
    write_aipc(My, AWD, qlSize);

    while(1){
/* Increments the PRIMARY semaphore immediately, i.e. lets it run
   and proceeds to the next line or exits if semop fails */
        if (semop(SemID, &buf0, 1) < 0) break;
        if (semop(SemID, &bufN, 1) < 0) break;

        AVARS = (struct A_vars *)ShmAd0;
        AARRS = (struct A_arrs *)ShmAd1;
        IOQL = (struct A_ql_IO *)ShmAdr;
// Fill IOQL with process information
	IOQL->My = My;

// Fill IOQL with Fortran-interface output arrays
/* Call Fortran function */
#ifdef qlk
        qlk_interf_(
#endif
#ifdef tglf
        tglf_interf_(
#endif
#ifdef neo
        neo_interf_(
#endif
           /* input */
	      &(IOQL->jrho_beg),
              &(IOQL->jrho_end),
              &(AVARS->na1),
              &(AVARS->na1n),
              &(AVARS->na1e),
              &(AVARS->na1i),
              &(AVARS->btor),
              &(AVARS->rtor),
              &(AVARS->amj),
              &(AVARS->zmj),
              &(AVARS->aim1),
              &(AVARS->aim2),
              &(AVARS->aim3),
              &(AARRS->ne),
              &(AARRS->te),
              &(AARRS->ni),
              &(AARRS->ndeut),
              &(AARRS->ntrit),
              &(AARRS->niz1),
              &(AARRS->niz2),
              &(AARRS->ti),
              &(AARRS->zef),
              &(AARRS->zim1),
              &(AARRS->amain),
              &(AARRS->mu),
              &(AARRS->rho),
              &(AARRS->ametr),
              &(AARRS->shif),
              &(AARRS->elon),
              &(AARRS->tria),
              &(AARRS->er),
              &(AARRS->nibm),
              &(AARRS->g11),
              &(AARRS->vpol),
              &(AARRS->vrs),
              &(AARRS->vtor),
              &(AARRS->shear),
              &(AARRS->pblon),
              &(AARRS->pbper),
              &(AARRS->pfast),
              &(AARRS->niz3),
              &(AARRS->zim2),
              &(AARRS->zim3),
              &(AARRS->zimpt),
              &(AARRS->nimpt),
              &(AARRS->aimpt),
/* output */
              &(IOQL->chi),
              &(IOQL->che),
              &(IOQL->dif),
              &(IOQL->vin),
              &(IOQL->dph),
              &(IOQL->dpl),
              &(IOQL->dpr),
              &(IOQL->xtb),
              &(IOQL->egm),
              &(IOQL->gam),
              &(IOQL->gm1),
              &(IOQL->gm2),
              &(IOQL->om1),
              &(IOQL->om2),
              &(IOQL->fr1)
          );

        swatch(&(My.CPUse));
/* If SemID exists then lock myself, otherwise, exit */
        Mysemun.val = 0;
        if (semctl(SemID, My.OrdNr, SETVAL, Mysemun) < 0) break;
    }

    if (errno != EIDRM && errno != EINVAL){
        printf(">>> %s >>> Unrecognised sem error: errno = %d\n",
            My.Path, errno);
        printf("EACCES = %d, EFAULT = %d, ERANGE = %d\n",
            EACCES, EFAULT, ERANGE);
        printf("EFBIG=%d, EINTR=%d, EAGAIN=%d, E2BIG=%d\n",
            EFBIG, EINTR, EAGAIN, E2BIG);
     }
    for(i=0; i < My.OrdNr; i++){
        printf("     ");
    }
    printf("Process # %d normal exit: ", My.OrdNr);
    swatch(&(My.CPUse));
    printf("CPUse %g\n", My.CPUse);
    exit(0);

}
