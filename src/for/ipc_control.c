#include "Astra.h"

int A_NA1;
int N_ARR_IN, N_ARR_OUT;

int read_aipc(int*);

char AWD[128];
char A_ipc_file[128];
key_t my_key;

int A_SemID = 0;
int A_Nsems = 0;        // #semaphores
void **A_ShmAdr = NULL;

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
void sbp2astra_(int* jrho_beg, int* jrho_end, int* jsbp, double* mem){
    int j, jarr, nrho;
    nrho = *jrho_end + 1 - *jrho_beg;
    double* sbp_out = (double *)((char *)A_ShmAdr[*jsbp-1]);
    for (j=*jrho_beg-1; j <= *jrho_end-1; j++){
        for (jarr=0; jarr<N_ARR_OUT; jarr++){
            mem[j+1+jarr*A_NA1] = sbp_out[j+jarr*nrho];
        }
    }
    return;
}

/*---------------------------------------------------
  Get PID and key for the Astra main process
  Create and initialize a set of A_Nsems semaphores
  Assign NA1 (= *Ngrid) to A_NA1
  Allocate two shared memory segments for Astra datasets
*/
int initialise_ipc_(int* Ngrid, int *n_sbp_arr_in, int *n_sbp_arr_out, int* Nsub, char* equ_file, char* exp_file){
  int var_size, arr_size, j, A_ShmID_avars, A_ShmID_aarrs;
    pid_t A_PID = 0;
    FILE *A_IPC;
    char hostname[128], ASTRA_task[128], A_equ_file[32], A_exp_file[32];
    time_t hold_time;
    static union semun Mysemun;

    A_Nsems = *Nsub + 1;
    A_ShmAdr = malloc(A_Nsems * sizeof(*A_ShmAdr));
    getcwd(AWD, sizeof(AWD));
    snprintf(A_equ_file, sizeof(A_equ_file), "%s", equ_file);
    snprintf(A_exp_file, sizeof(A_exp_file), "%s", exp_file);
    trim_right(A_equ_file);
    trim_right(A_exp_file);

// Collecting data
    A_PID = getpid();
// Define the absolute path name of Astra executable ASTRA_task
    snprintf(ASTRA_task, sizeof(ASTRA_task), "%s/bin/%s.exe", AWD, A_equ_file);
    my_key = ftok( ASTRA_task, (int)A_PID);    // Get System V IPC key
    printf("ASTRA_task %s\n", ASTRA_task);
    if (my_key == -1){
        fprintf(stderr, "ipc_control: not able to create Key from ProcID\n");
        fprintf(stderr, "Probably wrong ATASK name parsed from tmp/astra.nml\n");
        fprintf(stderr, "ATASK name: %s\n", ASTRA_task);
        exit(1);
    }

    if (A_Nsems != 0){
// Create a set of A_Nsems semaphores
        A_SemID = semget(my_key, A_Nsems, 0660|IPC_CREAT);
// Initialize all semaphores in the set A_SemID as {0, 0, ...}
        Mysemun.val = 0;
        for(j=0; j < A_Nsems; j++){
            semctl(A_SemID, j, SETVAL, Mysemun);
        }
    }
    
// Write file tmp/<exp><equ>.ipc
    snprintf(A_ipc_file, sizeof(A_ipc_file), "%s/tmp/%s%s.ipc", AWD, A_exp_file, A_equ_file);
    A_IPC = fopen(A_ipc_file, "w");
    if (!A_IPC){
        fprintf(stderr, "Cannot open Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(0);
    }
    fprintf(A_IPC, " Astra task:  \"%s\"\n", ASTRA_task);
    fprintf(A_IPC, " Astra files:  \"%s\",  \"%s\"\n", A_exp_file, A_equ_file);
    gethostname(hostname, (size_t)32);

    hold_time = time(NULL);
    fprintf(A_IPC, " Astra@%s started on:  %s", hostname, ctime(&hold_time));

    if (A_Nsems == 0){
        fclose(A_IPC);
        return 0;
    }
    fprintf(A_IPC, " Astra(main):  PID = %d,  SemID = %d\n", (int)A_PID, A_SemID);
    A_NA1 = *Ngrid;
    N_ARR_IN  = *n_sbp_arr_in;
    N_ARR_OUT = *n_sbp_arr_out;

    var_size = sizeof(struct A_vars);
    arr_size = A_NA1*N_ARR_IN*sizeof(double);
// Allocate shared memory segment for AVARS, AARRS
    A_ShmID_avars = shmget(my_key  , var_size, 0660|IPC_CREAT|IPC_EXCL);
    A_ShmID_aarrs = shmget(my_key+1, arr_size, 0660|IPC_CREAT|IPC_EXCL);
    A_ShmAdr_avars = shmat(A_ShmID_avars, NULL, 0);
    A_ShmAdr_aarrs = shmat(A_ShmID_aarrs, NULL, 0);
    fprintf(A_IPC, " Astra_inout: ShmID(vars):%12d%12d\n", A_ShmID_avars, var_size);
    fprintf(A_IPC, " Astra_inout: ShmID(arrs):%12d%12d\n", A_ShmID_aarrs, arr_size);
    fprintf(A_IPC, "        PID      ShmID       ShmSize     Process\n");
    fclose(A_IPC);
    return 0;
}

/*---------------------------------------------------------------------
  Set (lock) the primary semaphore to -(Number_of_processes)
  Launch parallel subprocesses
*/
int send_ipc_jobs_(int* Nsub, int *stringLen, char *subs, int* jbeg_arr, int* jend_arr){
    char jobString[400], path[32];
    int j, i;

    if (A_Nsems <= *Nsub){
        fprintf(stderr, " >>> ERROR >>> Incomplete semaphore set\n");
        return 1;
    }

    for (j=0; j<*Nsub; j++) {
        snprintf(path, sizeof(path), "%s", &subs[*stringLen * j]);
// Sending main (e.g. "tglfi"), only once per subprocess
        snprintf(jobString, sizeof(jobString), "%s/%s %s %d %d %d %d %d %d %d&",
		 AWD, path, A_ipc_file, my_key, j + 1, jbeg_arr[j], jend_arr[j], A_NA1, N_ARR_IN, N_ARR_OUT);
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
    if (read_aipc(Nsub)){
        fprintf(stderr, "Error in input data interpretation (function read_aipc)\n");
        exit(j);
    }
    return 0;
}

/*------------------------------------------------------*/
void WhatSem(){
    int j;
    if (A_Nsems == 0) return;
    ushort semarray[A_Nsems];
    union semun Mysemun;
    Mysemun.array = &semarray[0];
    semctl(A_SemID, 0, GETALL, Mysemun);

    printf(" Semaphore set = {");
    for (j=0; j<A_Nsems; j++) printf("%d, ", semarray[j]);
    printf("}\n");
    return;
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
                printf("\n");
                WhatSem();
                printf("ERANGE error:  sem_op = %d,   SEMVMX = ?\n\n", buf0.sem_op);
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
int read_aipc(int* Nsub){
    FILE *A_IPC;
    char ipcLine[128], childName[128];
    int i, j, ID, ShmID, ShmLen;
    char A_ChildName[A_Nsems][128]; // Child process name
    int A_ShmLen[A_Nsems];          // Child Shmem segment length
    int A_ShmID[A_Nsems];           // Child Shmem ID
    int A_ChildID[A_Nsems]; // Child process ID

    A_IPC = fopen(A_ipc_file, "r");
    if (!A_IPC){
        fprintf(stderr, "Cannot open Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(1);
    }
// Skip 7 header lines
    for (j=0; j<7; j++){
        fgets(ipcLine, sizeof(ipcLine), A_IPC);
    }
    i = 0;
    while (fscanf(A_IPC, "%12d%12d%12d%s", &ID, &ShmID, &ShmLen, childName) != EOF){
	if (i >= A_Nsems) {
            break;  // Too many subprocesses lines in IPC file
        }
        A_ChildID[i] = ID;
        A_ShmLen[i]  = ShmLen;
        A_ShmID[i]   = ShmID;
        A_ShmAdr[i]  = shmat(ShmID, NULL, 0);
        strcpy(&A_ChildName[i][0], childName);
        i++;
    }
    fclose(A_IPC);

    if ( i == *Nsub ) return 0;

    fprintf(stderr, " >>> File \"%s\" error >>> Missing SBP(%d)\n", A_ipc_file, i+1);
    fprintf(stderr, "\n  File \"%s\" contents:\n", A_ipc_file);
    fprintf(stderr, "\n  File processing:\n");
    fprintf(stderr, "A_Nsems = %d, A_ShmNum = %d,  Nsub = %d\n", A_Nsems, i-1, *Nsub);
    for (j=0; j <= *Nsub; j++){
        fprintf(stderr, "%12d%12d%12d%12d  %s\n", j, A_ChildID[j],
                A_ShmLen[j], A_ShmID[j], A_ChildName[j]);
    }
    return 1;
}

/*----------- Fill subprocess input scalars from ASTRA -------------*/
int setvars_(double* DEVAR, int* NBOUND){

    if (A_Nsems == 0){
        fprintf(stderr, " >>> setvars >>> Illegal call: Semaphores are not created\n");
        return 0;
    }

    AVARS = (struct A_vars *)A_ShmAdr_avars;
    AVARS->aim1 = *(DEVAR + 2);
    AVARS->aim2 = *(DEVAR + 3);
    AVARS->aim3 = *(DEVAR + 4);
    AVARS->amj  = *(DEVAR + 5);
    AVARS->btor = *(DEVAR + 7);
    AVARS->rtor = *(DEVAR + 27);
    AVARS->zmj  = *(DEVAR + 36);
    AVARS->na1n = *(NBOUND);
    AVARS->na1e = *(NBOUND + 1);
    AVARS->na1i = *(NBOUND + 2);
    return 0;
}

/*------------- Fill subprocess input arrays from ASTRA --------------*/
int fill_arr2shm_(double* sbp_in){

    double* sbp_input = (double *)((char *)A_ShmAdr_aarrs);
    size_t num_elements = A_NA1 * N_ARR_IN;
    memcpy(sbp_input, sbp_in, num_elements * sizeof(double));

    return 0;
}
