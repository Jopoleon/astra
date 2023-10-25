/*---------------------------------------------------------------------*/
#define C_MAIN 0
/*------------------------------------------------------------------*/
/* Note! It is expected that the variable C_MAIN is defined in 
         the C-file where this file __FILE__ is included            */
/* The file includes functions (double)`swatch' and 
   the starting part of the function (int)`sbp2shm' that
   acquires information about AWD, process wd, id, name, key, 
   ord. number, allocated shared memory ID and semaphore ID.        */
/*------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <errno.h>
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

#include "A_proc.h"

void a_stop_();

extern INT_  A_NB1;
int   SemID ,  ShMid0,  ShMid1;
void *ShmAd0, *ShmAd1, *ShmAdr;
void qlk_interf_();
void neo_interf_();
void tglf_interf_();
void SP_stamp();
int read_aipc();
int write_aipc();

#define NC1 A_NB1

/*---------------------------------------------------------------------*/
#if C_MAIN != 0
int main(argc, argv, envir)   /* C main __FILE__  is used */
int argc;
char *argv[], *envir[];
#else                         /* Fortran main is used */
   /* Note description of arg2, arg3, arg4 should be compatible with
      INT_ or int
   */

int sbp2shm_(char* arg0, char* arg1, int* arg2, int* arg3, int* arg4)
#endif
{
    int i, lS;
    static union semun Mysemun;
// sembuf members: {sem_num,sem_op,sem_flag};
    static struct sembuf buf0 = {0, 1, IPC_NOWAIT};
    static struct sembuf bufN = {1,-1, ~SEM_UNDO&~IPC_NOWAIT};
    static char whoami[32];
    static char AWD[96];
#include "A_vars.h"
#include "A_proc.h"
    static struct A_proc_info Mama, My;

    swatch(&(My.CPUse));
/* Analyze the calling command string. Get own PID and name. */
    My.Pid = getpid();
    getcwd(My.Path, (size_t)64);
    strcat(My.Path, "/");

#if C_MAIN != 0
    strcat(My.Path,argv[0]+2);     /* C main __FILE__  is used */
    if (strrchr(argv[0],'/') != NULL){
        sscanf(strrchr(argv[0], '/')+1,"%s",whoami);
    }
    else{
        sscanf(argv[0], "%s", whoami);
    }
    sscanf(argv[1], "%s", Mama.Path);
    ++argv;
    ++argv;
    Mama.Pid = (pid_t)atoi(*argv);
    ++argv;
    Mama.Key  = (key_t)atoi(*argv);
    ++argv;
    My.OrdNr = atoi(*argv);
    printf("C main: %s\n", My.Path);
#else            /* Fortran main is used */
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
#endif

/* Associate My semaphore with the ordinal process number */
    bufN.sem_num = My.OrdNr;
/* Get semaphore and shmem IDs. */
    SemID  = semget(Mama.Key, 0, 0660);
    ShMid0 = shmget((key_t)(Mama.Key+0), 0, 0660);
    ShMid1 = shmget((key_t)(Mama.Key+1), 0, 0660);
    ShmAd0 = shmat(ShMid0, NULL, 0);
    AVARS = (struct A_vars *)ShmAd0;
    A_NB1 = AVARS->nb1;
    ShmAd1 = shmat(ShMid1, NULL, 0);
    SP_stamp(Mama, My, SemID);
    strcpy(AWD, Mama.Path);
    if (strstr(AWD, "bin/") == NULL){
        printf("SBP launch string error\n");
        a_stop_();
    }
    else{
        *strstr(AWD, "bin/") = '\0';
    }
    My.Key = ftok( My.Path, (int)My.Pid);

#include "A_arrs.h"
#include "A_ql_IO.h"

    lS = sizeof(struct A_ql_IO);

/*------------------------------------
  Create My shared memory segment
*/
    My.ShMid = shmget(My.Key, lS, 0660|IPC_CREAT);
/* Attach My shared memory to the process
   My shared memory segment starts at ShmAdr */
    ShmAdr = shmat(My.ShMid, NULL, 0);
    i = My.OrdNr;
    i = write_aipc(My, AWD, &lS);

    while(1){
/* Increments the PRIMARY semaphore immediately, i.e. lets it run
   and proceeds to the next line or exits if semop fails */
        if (semop(SemID, &buf0, 1) < 0) break;
        if (semop(SemID, &bufN, 1) < 0) break;

        AVARS = (struct A_vars *)ShmAd0;
        AARRS = (struct A_arrs *)ShmAd1;
        IOQL = (struct A_ql_IO *)ShmAdr;
        IOQL->My = My;

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
              &(IOQL->is), 
              &(IOQL->ie), 
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
              &((*IOQL).chi), 
              &((*IOQL).che), 
              &((*IOQL).dif), 
              &((*IOQL).vin), 
              &((*IOQL).dph), 
              &((*IOQL).dpl), 
              &((*IOQL).dpr), 
              &((*IOQL).xtb), 
              &((*IOQL).egm), 
              &((*IOQL).gam), 
              &((*IOQL).gm1), 
              &((*IOQL).gm2), 
              &((*IOQL).om1), 
              &((*IOQL).om2), 
              &((*IOQL).fr1)
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
