#include "Astra.h"

void qlk_interf_(int*, int*, int*, int*, int*, int*, int*, int*,
		 double*, double*, double*, double*, double*, double*, double*,
		 double*, double*);
void neo_interf_(int*, int*, int*, int*, int*, int*, int*, int*,
		 double*, double*, double*, double*, double*, double*, double*,
		 double*, double*);
void tglf_interf_(int*, int*, int*, int*, int*, int*, int*, int*,
		 double*, double*, double*, double*, double*, double*, double*,
		 double*, double*);

/*---------------------------------------------------------------------*/
/* acquires information about AWD, process wd, id, name, key,
   ord. number, allocated shared memory ID and semaphore ID.
*/

int main(int argc, char *argv[]) {

    extern int A_NA1;
    void *ShmAd0, *ShmAd1, *ShmAdr;
    int i, qlSize, N_ARR_IN, N_ARR_OUT, jrho_beg, jrho_end;
    int SemID, ShmId0, ShmId1;
    int ProcOrdNr, ProcShmId;
    key_t Key_in, ProcKey;
    pid_t ProcPid;
    char ProcPath[128];
    static union semun Mysemun;
// sembuf members: {sem_num,sem_op,sem_flag};
    static struct sembuf buf0 = {0, 1, IPC_NOWAIT};
    static struct sembuf bufN = {1,-1, ~SEM_UNDO&~IPC_NOWAIT};
    static char A_ipc_file[128];
    FILE *A_IPC;

/* Analyze the calling command string. Get own PID and name. */
    ProcPid = getpid();
    strcpy(ProcPath, argv[0]);

    sscanf(argv[1], "%s", A_ipc_file);
    Key_in = (key_t)atoi(argv[2]);
    ProcOrdNr = atoi(argv[3]);
    jrho_beg  = atoi(argv[4]);
    jrho_end  = atoi(argv[5]);
    A_NA1     = atoi(argv[6]);
    N_ARR_IN  = atoi(argv[7]);
    N_ARR_OUT = atoi(argv[8]);
    printf("Fortran main: %s %d %3d %3d %d %d %d\n", ProcPath, ProcPid, jrho_beg, jrho_end, A_NA1, N_ARR_IN, N_ARR_OUT);

/* Associate My semaphore with the ordinal process number */
    bufN.sem_num = ProcOrdNr;
/* Get semaphore and shmem IDs. */
    SemID  = semget(Key_in, 0, 0660);
    ShmId0 = shmget((key_t)(Key_in+0), 0, 0660);
    ShmId1 = shmget((key_t)(Key_in+1), 0, 0660);
    ShmAd0 = shmat(ShmId0, NULL, 0);
    ShmAd1 = shmat(ShmId1, NULL, 0);
    ProcKey = ftok( ProcPath, (int)ProcPid);
    qlSize = (jrho_end + 1 - jrho_beg)*N_ARR_OUT*sizeof(double);

/*------------------------------------
  Create My shared memory segment
*/
    ProcShmId = shmget(ProcKey, qlSize, 0660|IPC_CREAT);
/* Attach My shared memory to the process
   My shared memory segment starts at ShmAdr */
    ShmAdr = shmat(ProcShmId, NULL, 0);

// Append process info to A_ipc_file
    A_IPC = fopen(A_ipc_file, "a");
    if (!A_IPC){
        printf("Cannot open existing Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(0);
    }
    fprintf(A_IPC, "%12d%12d%12d   %s\n", ProcPid, ProcShmId, qlSize, ProcPath);
    fclose(A_IPC);

    while(1){
/* Increments the PRIMARY semaphore immediately, i.e. lets it run
   and proceeds to the next line or exits if semop fails */
        if (semop(SemID, &buf0, 1) < 0) break;
        if (semop(SemID, &bufN, 1) < 0) break;

        AVARS = (struct A_vars *)ShmAd0;
	double* sbp_in = (double *)((char *)ShmAd1);
        double* sbp_out = (double *)((char *)ShmAdr); // IPC subprocess output

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
            &jrho_beg,
            &jrho_end,
            &N_ARR_IN,
            &N_ARR_OUT,
            &A_NA1,
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
	    sbp_in,
/* output */
            sbp_out
          );

/* If SemID exists then lock myself, otherwise, exit */
        Mysemun.val = 0;
        if (semctl(SemID, ProcOrdNr, SETVAL, Mysemun) < 0) break;
    }

    printf("Process # %d normal exit:\n", ProcOrdNr);
    exit(0);
}
