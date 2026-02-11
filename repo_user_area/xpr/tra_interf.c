#include "Astra.h"

void qlk_interf_( int*, int*, double*, double*, double*);
void neo_interf_( int*, int*, double*, double*, double*);
void tglf_interf_(int*, int*, double*, double*, double*);

/*---------------------------------------------------------------------*/
/* acquires information about AWD, process wd, id, name, key,
   ord. number, allocated shared memory ID and semaphore ID.
*/

int main(int argc, char *argv[]) {

    void *ShmAd0, *ShmAd1, *ShmAd2, *ShmAdr;
    int j, outSize, N_ARR_OUT, J_PROC, N_CHUNK;
    int SemID, ShmId0, ShmId1, ShmId2;
    int ProcShmId;
    key_t Key_in, ProcKey;
    pid_t ProcPid;
    char ProcPath[128];
    static union semun Mysemun;
// sembuf members: {sem_num,sem_op,sem_flag};
    static struct sembuf buf0 = {0, 1, IPC_NOWAIT};
    static struct sembuf bufN = {1,-1, ~SEM_UNDO&~IPC_NOWAIT};
    static char ipc_file[128];
    FILE *IPCa;
    FILE *IPCr;

/* Analyze the calling command string. Get own PID and name. */
    ProcPid = getpid();
    strcpy(ProcPath, argv[0]);

    sscanf(argv[1], "%s", ipc_file);
    Key_in = (key_t)atoi(argv[2]);
    J_PROC    = atoi(argv[3]);
    N_CHUNK   = atoi(argv[4]);
    N_ARR_OUT = atoi(argv[5]);
    ShmId0 = atoi(argv[6]);
    ShmId1 = atoi(argv[7]);
    ShmId2 = atoi(argv[8]);
    printf("Fortran main: %s %d %3d %3d %d\n", ProcPath, ProcPid, J_PROC, N_CHUNK, N_ARR_OUT);

/* Associate My semaphore with the ordinal process number */
    bufN.sem_num = J_PROC;
/* Get semaphore and shmem IDs. */
    SemID  = semget(Key_in, 0, 0660);
    ShmAd0 = shmat(ShmId0, NULL, 0);
    ShmAd1 = shmat(ShmId1, NULL, 0);
    ShmAd2 = shmat(ShmId2, NULL, 0);
    ProcKey = ftok(ProcPath, (int)ProcPid);
    outSize = N_CHUNK*N_ARR_OUT*sizeof(double);

/*------------------------------------
  Create My shared memory segment
*/
    ProcShmId = shmget(ProcKey, outSize, 0660|IPC_CREAT);
/* Attach My shared memory to the process
   My shared memory segment starts at ShmAdr */
    ShmAdr = shmat(ProcShmId, NULL, 0);

// Append process info to ipc_file
    IPCa = fopen(ipc_file, "a");
    if (!IPCa){
        printf("Cannot open existing Astra IPC file: \"%s\"\n", ipc_file);
        exit(0);
    }
    fprintf(IPCa, "%10d %10d %10d\n", J_PROC, ProcPid, ProcShmId);
    fclose(IPCa);

    int* dim_in = (int *)((char *)ShmAd0); // constant at all time steps

    while(1){
/* Increments the PRIMARY semaphore immediately, i.e. lets it run
   and proceeds to the next line or exits if semop fails */
        if (semop(SemID, &buf0, 1) < 0) break;
        if (semop(SemID, &bufN, 1) < 0) break;
        double* scal_in  = (double *)((char *)ShmAd1);
	double* prof_in  = (double *)((char *)ShmAd2);
        double* prof_out = (double *)((char *)ShmAdr); // IPC subprocess output

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
       	    &J_PROC, dim_in, scal_in, prof_in, // input
            prof_out // output
        );

/* If SemID exists then lock myself, otherwise, exit */
        Mysemun.val = 0;
        if (semctl(SemID, J_PROC, SETVAL, Mysemun) < 0) break;
    }

    printf("Process # %d normal exit:\n", J_PROC);
    exit(0);
}
