#include "Astra.h"

int read_aipc(int*, int*, char*);

int A_SemID = 0;
void **A_ShmAdr = NULL;
void *A_ShmAdr_adims;
void *A_ShmAdr_avars;
void *A_ShmAdr_aarrs;

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
void sbp2astra_(int* jsbp, int *nchunk, int *n_sbp_arr_out, double* mem){

    int j, jproc, jarr, N_CHUNK, N_ARR_OUT;
    jproc = *jsbp - 1;
    N_CHUNK = *nchunk;
    N_ARR_OUT = *n_sbp_arr_out;
    double* prof_out = (double *)((char *)A_ShmAdr[jproc]);

    for (j=0; j<N_CHUNK; j++){
        for (jarr=0; jarr<N_ARR_OUT; jarr++){
            mem[jarr + (j + jproc*N_CHUNK) * N_ARR_OUT] = prof_out[jarr + j*N_ARR_OUT];
        }
    }
    return;
}

/*---------------------------------------------------
  Get PID and key for the Astra main process
  Create and initialize a set of A_Nsems semaphores
  Assign NA1 (= *Ngrid) to N_RHO
  Allocate two shared memory segments for Astra datasets
*/
int initialise_ipc_(int* Ngrid, int* Ndims, int* Nscalars, int *n_sbp_arr_in, int *n_sbp_arr_out, int *nchunk, int* Nsub, int* offset, char *subName, char* ipc_file, char* astra_task){
  
    int N_DIMS, N_SCALARS, N_RHO, N_ARR_IN, N_ARR_OUT, N_CHUNK, dim_size, var_size, arr_size, A_Nsems;
    int i, j, A_ShmID_adims, A_ShmID_avars, A_ShmID_aarrs;
    key_t my_key;
    pid_t A_PID=0;
    FILE *A_IPCw;
    char hostname[128], A_astra_task[128], A_ipc_file[128], A_sub_name[64];
    time_t hold_time;
    static union semun Mysemun;

    A_Nsems = *Nsub + 1;
    N_DIMS = *Ndims;
    N_SCALARS = *Nscalars;
    N_RHO = *Ngrid;
    N_ARR_IN  = *n_sbp_arr_in;
    N_ARR_OUT = *n_sbp_arr_out;
    N_CHUNK = *nchunk;

    A_ShmAdr = malloc(*Nsub * sizeof(*A_ShmAdr));
    snprintf(A_sub_name, sizeof(A_sub_name), "%s", subName);
    snprintf(A_ipc_file, sizeof(A_ipc_file), "%s", ipc_file);
    snprintf(A_astra_task, sizeof(A_astra_task), "%s", astra_task);
    trim_right(A_ipc_file);
    trim_right(A_astra_task);
    trim_right(A_sub_name);

// Collecting data
    A_PID = getpid();
// Define the absolute path name of Astra executable ASTRA_task
    my_key = ftok( A_astra_task, (int)A_PID);    // Get System V IPC key
    printf("ASTRA_task %s %d %d\n", A_astra_task, my_key, (int)A_PID);
    if (my_key == -1){
        fprintf(stderr, "ipc_control: not able to create Key from ProcID\n");
        fprintf(stderr, "Probably wrong ATASK name parsed from tmp/astra.nml\n");
        fprintf(stderr, "ATASK name: %s\n", A_astra_task);
        exit(1);
    }

// Create a set of A_Nsems semaphores
    A_SemID = semget(my_key, A_Nsems, 0660|IPC_CREAT);
// Initialize all semaphores in the set A_SemID as {0, 0, ...}
    Mysemun.val = 0;
    for(j=0; j<A_Nsems; j++){
        semctl(A_SemID, j, SETVAL, Mysemun);
    }

    dim_size = N_DIMS*sizeof(int);
    var_size = N_SCALARS*sizeof(double);
    arr_size = N_RHO*N_ARR_IN*sizeof(double);
// Allocate shared memory segment for AVARS, AARRS
    A_ShmID_adims = shmget(my_key  , dim_size, 0660|IPC_CREAT|IPC_EXCL);
    A_ShmID_avars = shmget(my_key+1, var_size, 0660|IPC_CREAT|IPC_EXCL);
    A_ShmID_aarrs = shmget(my_key+2, arr_size, 0660|IPC_CREAT|IPC_EXCL);
    A_ShmAdr_adims = shmat(A_ShmID_adims, NULL, 0);
    A_ShmAdr_avars = shmat(A_ShmID_avars, NULL, 0);
    A_ShmAdr_aarrs = shmat(A_ShmID_aarrs, NULL, 0);

// Write file tmp/<exp><equ>.ipc
    A_IPCw = fopen(A_ipc_file, "w");
    if (!A_IPCw){
        fprintf(stderr, "Cannot open Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(0);
    }
    fprintf(A_IPCw, "Astra task:  \"%s\"\n", A_astra_task);
    gethostname(hostname, (size_t)32);

    hold_time = time(NULL);
    fprintf(A_IPCw, "Astra@%s started on:  %s", hostname, ctime(&hold_time));
    fprintf(A_IPCw, "Astra(main):\n");
    fprintf(A_IPCw, "  PID   = %d\n", (int)A_PID);
    fprintf(A_IPCw, "  SemID = %d\n", A_SemID);
    fprintf(A_IPCw, "  Key   = %d\n", (int)my_key);
    fprintf(A_IPCw, "Dims, ShmId, size:%12d%12d\n", A_ShmID_adims, dim_size);
    fprintf(A_IPCw, "Vars: ShmID, size:%12d%12d\n", A_ShmID_avars, var_size);
    fprintf(A_IPCw, "Arrs: ShmID, size:%12d%12d\n", A_ShmID_aarrs, arr_size);
    fprintf(A_IPCw, "         PID       ShmID\n");
    fclose(A_IPCw);

/*---------------------------------------------------------------------
  Set (lock) the primary semaphore to -(Number_of_processes)
  Launch parallel subprocesses
*/

    char jobString[400], AWD[128];

    getcwd(AWD, sizeof(AWD));
    for (j=0; j<*Nsub; j++) {
// Sending main (e.g. "tglfi"), only once per subprocess
        snprintf(jobString, sizeof(jobString), "%s/%s %s %d %d %d %d &",
		 AWD, A_sub_name, A_ipc_file, my_key, j + 1, N_CHUNK, N_ARR_OUT);
        i = system(jobString);

// Wait until child increments semaphore 0
        struct sembuf bufj = {0, -1, ~IPC_NOWAIT};
        if (semop(A_SemID, &bufj, 1) == -1) {
            perror("semop failed");
            return j + 1;
        }

        if (i == -1) return j + 1;
    }

// Read process Shm addresses and size from A_ipc_file
    if (read_aipc(Nsub, offset, ipc_file)){
        fprintf(stderr, "Error in input data interpretation (function read_aipc)\n");
        exit(j);
	}
    return 0;
}

/*--------------------- Unlock subprocess ----------------------------*/
void unlock_sbp_(int* jsbp){
    auto struct sembuf bufJ = {*jsbp, 1, IPC_NOWAIT};
    --buf0.sem_op;   // Each call decrements Sem0 value by 1
// Increment semval # sem_num=*n, Open subprocess 
    semop(A_SemID, &bufJ, 1);
}

/*---------------------------------------------------------------------*/
/* Compare Sem0 value with buf0.sem_op ( == -A_Nsems) and
 wait until all subprocesses increment Sem0 by 1
        so that Sem0 reaches value A_Nsems
*/
int wait4all_(){
    static struct timespec timeout = {0, 100000000};   // timeout = .1 sec

/*
  Note! Each call increments value of semadj by 1.
  Overflow occurs when the number of successive calls exceeds 32767
*/
    while(1){
        if (semtimedop(A_SemID, &buf0, 1, &timeout) == 0) {
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

/*------------------------------------------------------------
  Read IPC file after launching all processes and
    (1) fill arrays of Child_process_IDs, Child_shmem_lengths/IDs,
    (2) attach child process shmem segments to the main process memory
*/
int read_aipc(int* Nsub, int *offset, char *ipc_file){

    char A_ipc_file[128];
    snprintf(A_ipc_file, sizeof(A_ipc_file), "%s", ipc_file);
    trim_right(A_ipc_file);
    FILE *A_IPCr = fopen(A_ipc_file, "r");
    int i, j, c, ID, ShmID;
    
    if (!A_IPCr) {
        perror(A_ipc_file);
        exit(EXIT_FAILURE);
    }

/* Skip header lines */
    for (j=0; j < *offset; j++) {
        while ((c = fgetc(A_IPCr)) != '\n' && c != EOF);
    }

/* Read Proc IDs and assign address */
    i = 0;
    while (i < *Nsub &&
        fscanf(A_IPCr, "%12d%12d", &ID, &ShmID) == 2) {
        A_ShmAdr[i] = shmat(ShmID, NULL, 0);
	i++;
    }
    fclose(A_IPCr);
    return 0;
}

/*----------- Fill subprocess input dims from ASTRA -------------*/
int fill_dim2shm_(int* Ndims, int* dims_in){

    int* dims_input = (int *)((char *)A_ShmAdr_adims);
    int N_DIMS = *Ndims;
    int dim_size = N_DIMS * sizeof(int);
    memcpy(dims_input, dims_in, dim_size);

    return 0;
}

/*----------- Fill subprocess input scalars from ASTRA -------------*/
int fill_var2shm_(int* Nscalars, double* scal_in){

    double* scal_input = (double *)((char *)A_ShmAdr_avars);
    int N_SCALARS = *Nscalars;
    int var_size = N_SCALARS*sizeof(double);
    memcpy(scal_input, scal_in, var_size);

    return 0;
}

/*------------- Fill subprocess input arrays from ASTRA --------------*/
int fill_arr2shm_(int* Ngrid, int* n_sbp_arr_in, double* prof_in){

    double* prof_input = (double *)((char *)A_ShmAdr_aarrs);
    int N_RHO = *Ngrid;
    int N_ARR_IN = *n_sbp_arr_in;
    size_t num_elements = N_RHO * N_ARR_IN;
    memcpy(prof_input, prof_in, num_elements * sizeof(double));

    return 0;
}
