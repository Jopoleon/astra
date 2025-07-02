#include "Astra.h"

extern INT_ A_NA1;
INT_ NQL;

void a_stop_();
int   SemID ,  ShMid0,  ShMid1;
void *ShmAd0, *ShmAd1, *ShmAdr;
double swatch_();
void qlk_interf_();
void neo_interf_();
void tglf_interf_();
void write_aipc();

/*---------------------------------------------------------------------*/
/* acquires information about AWD, process wd, id, name, key,
   ord. number, allocated shared memory ID and semaphore ID.
   Note description of arg2, arg3, arg4 should be compatible with INT_ or int
*/

int sbp2shm_(char* arg0, char* arg1, int* arg2, int* arg3, int* arg4)
{
    int i, qlSize;
    double watch;
    static union semun Mysemun;
// sembuf members: {sem_num,sem_op,sem_flag};
    static struct sembuf buf0 = {0, 1, IPC_NOWAIT};
    static struct sembuf bufN = {1,-1, ~SEM_UNDO&~IPC_NOWAIT};
    static char whoami[32];
    static char AWD[96];
    static struct A_proc_info Mama, My;

    swatch_(&(My.CPUse));
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
    NQL   = AVARS->n_ql;
#include "A_arrs.h"
    strcpy(AWD, Mama.Path);
    if (strstr(AWD, "bin/") == NULL){
        printf("SBP launch string error\n");
        a_stop_();
    }
    else{
        *strstr(AWD, "bin/") = '\0';
    }
    My.Key = ftok( My.Path, (int)My.Pid);
    qlSize = sizeof(struct A_ql_io) - sizeof(double) + NQL*A_NA1*sizeof(double);

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
        ql_io = (struct A_ql_io *)ShmAdr;
// Fill ql_io with process information
	ql_io->My = My;

// Fill ql_io with Fortran-interface output arrays
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
	      &(ql_io->jrho_beg),
              &(ql_io->jrho_end),
              &(AVARS->n_ql),
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
              &(ql_io->QLarrays)
          );

        swatch_(&(My.CPUse));
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
    swatch_(&(My.CPUse));
    printf("CPUse %g\n", My.CPUse);
    exit(0);

}
