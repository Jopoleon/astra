#include "Astra.h"

struct sembuf buf0 = {0, 0, ~SEM_UNDO&~IPC_NOWAIT};

/*---------------------------------------------------
  Trims a string
*/
void trim_right(char *str) {
    int i = strlen(str) - 1;
    while (i >= 0 && (str[i] == ' ' || str[i] == '\n' || str[i] == '\r' || str[i] == '\t'))
        str[i--] = '\0';
}

/*--------------------------------------------------------------------
  Reads the shared memory segment and stores the subprocess' output to an ASTRA fortran array
*/
void sbp2astra_(int* jsbp, int *nchunk, int *n_sbp_arr_out, char* ipcFile, double* mem){

    int j, J_PROC, jproc, ID, ShmID, jarr, N_CHUNK, N_ARR_OUT;
    char ipc_file[128];
    FILE *IPCr;
    char line[256];
 
    snprintf(  ipc_file, sizeof(ipc_file)  , "%s", ipcFile);
    trim_right(ipc_file);

    J_PROC = *jsbp - 1;
    N_CHUNK = *nchunk;
    N_ARR_OUT = *n_sbp_arr_out;
    void *ShmAdr = NULL;

    IPCr = fopen(ipc_file, "r");
    while (fgets(line, sizeof line, IPCr)) {
        char *p = line;
/* skip leading blanks and tabs */
        while (*p == ' ' || *p == '\t')
            p++;
/* skip line if first non-blank is not digit or sign */
        if (!isdigit((unsigned char)*p) && *p != '+' && *p != '-')
            continue;
/* parse integers */
        if (sscanf(p, "%d%d%d", &jproc, &ID, &ShmID) != 3)
            continue;
	if (jproc == J_PROC) break;
    }
    fclose(IPCr);
    
    ShmAdr = shmat(ShmID, NULL, 0);
    double* prof_out = (double *)((char *)ShmAdr);

    for (j=0; j<N_CHUNK; j++){
        for (jarr=0; jarr<N_ARR_OUT; jarr++){
            mem[jarr + (j + J_PROC*N_CHUNK) * N_ARR_OUT] = prof_out[jarr + j*N_ARR_OUT];
        }
    }
    return;
}

/*---------------------------------------------------
  Get PID and key for the Astra main process
  Create and initialize a set of Nsems semaphores
  Assign NA1 (= *Ngrid) to N_RHO
  Allocate two shared memory segments for Astra datasets
*/
int initialise_ipc_(int* Ngrid, int* Ndims, int* Nscalars, int *n_sbp_arr_in, int *n_sbp_arr_out, int* Nsub, char *subName, char* ipcFile, char* astraTask, int *SemID, int *ShmID_dims, int *ShmID_vars, int *ShmID_arrs, int *ipcId){

    int N_SUB, N_DIMS, N_SCALARS, N_RHO, N_ARR_IN, N_ARR_OUT, N_CHUNK, dim_size, var_size, arr_size, Nsems;
    int i, j, ID, ShmID;
    key_t my_key, key2;
    pid_t PID=0;
    FILE *IPCw;
    char hostname[128], AWD[128], astra_task[128], ipc_file[128], sub_name[64], jobString[400];
    time_t hold_time;
    static union semun Mysemun;

    N_SUB = *Nsub;
    Nsems = N_SUB + 1;
    N_DIMS = *Ndims;
    N_SCALARS = *Nscalars;
    N_RHO = *Ngrid;
    N_ARR_IN  = *n_sbp_arr_in;
    N_ARR_OUT = *n_sbp_arr_out;
    N_CHUNK = N_RHO/N_SUB;

    snprintf(  sub_name, sizeof(sub_name)  , "%s", subName);
    snprintf(  ipc_file, sizeof(ipc_file)  , "%s", ipcFile);
    snprintf(astra_task, sizeof(astra_task), "%s", astraTask);
    trim_right(ipc_file);
    trim_right(astra_task);
    trim_right(sub_name);
 
    PID = getpid();
// Define the absolute path name of Astra executable ASTRtask
    my_key = ftok( astra_task, (int)PID);    // Get System V IPC key
    dim_size = N_DIMS*sizeof(int);
    var_size = N_SCALARS*sizeof(double);
    arr_size = N_RHO*N_ARR_IN*sizeof(double);

// Collecting data
    if (my_key == -1){
        fprintf(stderr, "ipc_control: not able to create Key from ProcID\n");
        fprintf(stderr, "Probably wrong ATASK name parsed from tmp/astra.nml\n");
        fprintf(stderr, "ATASK name: %s\n", astra_task);
        exit(1);
    }

    key2 = my_key + (*ipcId)*3;
    printf("ASTRtask %s %d %d %d\n", astra_task, key2, (int)PID, dim_size);
// Create a set of Nsems semaphores
    *SemID = semget(key2, Nsems, 0660|IPC_CREAT);
// Initialize all semaphores in the set SemID as {0, 0, ...}
    Mysemun.val = 0;
    for(j=0; j<Nsems; j++){
        semctl(*SemID, j, SETVAL, Mysemun);
    }

// Allocate shared memory segment for VARS, ARRS
    *ShmID_dims = shmget(key2  , dim_size, 0660|IPC_CREAT|IPC_EXCL);
    *ShmID_vars = shmget(key2+1, var_size, 0660|IPC_CREAT|IPC_EXCL);
    *ShmID_arrs = shmget(key2+2, arr_size, 0660|IPC_CREAT|IPC_EXCL);

// Write file tmp/<exp><equ>.ipc
    IPCw = fopen(ipc_file, "w");
    if (!IPCw){
        fprintf(stderr, "Cannot open Astra IPC file: \"%s\"\n", ipc_file);
        exit(0);
    }
    fprintf(IPCw, "Astra task:  \"%s\"\n", astra_task);
    gethostname(hostname, (size_t)32);
    hold_time = time(NULL);
    fprintf(IPCw, "Astra@%s started on:  %s", hostname, ctime(&hold_time));
    fprintf(IPCw, "Astra(main):\n");
    fprintf(IPCw, "  PID   : %d\n", (int)PID);
    fprintf(IPCw, "  SemID : %d\n", *SemID);
    fprintf(IPCw, "  Key   : %d\n", (int)key2);
    fprintf(IPCw, "Dims | ShmId, size:%12d%12d\n", *ShmID_dims, dim_size);
    fprintf(IPCw, "Vars | ShmID, size:%12d%12d\n", *ShmID_vars, var_size);
    fprintf(IPCw, "Arrs | ShmID, size:%12d%12d\n", *ShmID_arrs, arr_size);
    fprintf(IPCw, "      ProcID       ShmID\n");
    fclose(IPCw);

/*---------------------------------------------------------------------
  Set (lock) the primary semaphore to -(Number_of_processes)
  Launch parallel subprocesses
*/

    getcwd(AWD, sizeof(AWD));
    for (j=0; j<*Nsub; j++) {
// Sending main (e.g. "tglfi"), only once per subprocess
        snprintf(jobString, sizeof(jobString), "%s/%s %s %d %d %d %d %d %d %d&",
		 AWD, sub_name, ipc_file, key2, j+1, N_CHUNK, N_ARR_OUT,
		 *ShmID_dims, *ShmID_vars, *ShmID_arrs);
        i = system(jobString);

/* Wait until child increments semaphore 0 */
        struct sembuf bufj = {0, -1, ~IPC_NOWAIT};
        if (semop(*SemID, &bufj, 1) == -1) {
            perror("semop failed");
            return j + 1;
        }
        if (i == -1) return j + 1;
    }

    return 0;
}

/*--------------------- Unlock subprocess ----------------------------*/
void unlock_sbp_(int* jsbp, int *SemID){
    auto struct sembuf bufJ = {*jsbp, 1, IPC_NOWAIT};
    --buf0.sem_op;   // Each call decrements Sem0 value by 1
// Increment semval # sem_num=*n, Open subprocess 
    semop(*SemID, &bufJ, 1);
}

/*---------------------------------------------------------------------*/
/* Compare Sem0 value with buf0.sem_op ( == -Nsems) and
 wait until all subprocesses increment Sem0 by 1
        so that Sem0 reaches value Nsems
*/
int wait4all_(int *SemID){
    static struct timespec timeout = {0, 100000000};   // timeout = .1 sec

/*
  Note! Each call increments value of semadj by 1.
  Overflow occurs when the number of successive calls exceeds 32767
*/
    while(1){
        if (semtimedop(*SemID, &buf0, 1, &timeout) == 0) {
            buf0.sem_op = 0;
            return 0;
        }
        switch (errno) {
            case EAGAIN: // Timeout: retry
                continue;
            case EIDRM:
                printf("The semaphore set was removed\n");
                break;
            case EINTR:
                printf("While blocked in this call, the process caught a signal\n");
                break;
            case ERANGE:
                printf("ERANGE:\t");
                break;
            case EINVAL:
                printf("EINVAL:\t");
                break;
            case EFAULT:
                printf("EFAULT:\t");
                break;
            case E2BIG:
                printf("E2BIG:\t");
                break;
            case EACCES:
                printf("EACCES:\t");
                break;
            case EFBIG:
                printf("EFBIG:\t");
                break;
            case ENOMEM:
                printf("ENOMEM:\t");
                break;
            default:
                printf(">>> PRIMARY >>> Unrecognised semtimedop error\n");
                break;
        }
        fprintf(stderr, ">>> PRIMARY >>> semtimedop error = %d\n", errno);
        break; 
    }
    return 1;
}

/*----------- Fill shared memory segment with integer array (dims) -------------*/
int fill_int_shm_(int* Nsize, int* int_in, int* ShmID){

    void *ShmAdr = shmat(*ShmID, NULL, 0);
    int* int_input = (int *)((char *)ShmAdr);
    int N_SIZE = *Nsize;
    int buf_size = N_SIZE*sizeof(int);
    memcpy(int_input, int_in, buf_size);

    return 0;
}

/*----------- Fill shared memory segment with double array (scal, prof) -------------*/
int fill_dbl_shm_(int* Nsize, double* dbl_in, int* ShmID){

    void *ShmAdr = shmat(*ShmID, NULL, 0);
    double* dbl_input = (double *)((char *)ShmAdr);
    int N_SIZE = *Nsize;
    int buf_size = N_SIZE*sizeof(double);
    memcpy(dbl_input, dbl_in, buf_size);

    return 0;
}
