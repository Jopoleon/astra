#include <unistd.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <signal.h>
#include <ctype.h>
#include <stddef.h>
#include "Astra.h"

INT_ A_NA1;
INT_ N_QL;

int semtimedop();
int read_aipc(INT_*, INT_*, char*);
void freeshm();
void to_tra_(INT_*, INT_*, INT_*);
void ot_tra_(INT_*, INT_*, INT_*, double*, double*);
void AllocateShmem(int);

char *AWD, *equmod, *DATA;
char A_ipc_file[132];
char ASTRA_task[132];
const char *A_log_file = "./tmp/astra.nml";
key_t my_key;

pid_t A_PID = 0;
int A_SemID = 0;
int A_Nsems = 0;       /* the number of semaphores */
int A_ShmNum = -1;
#define A_ShmShift 2
#define A_Nsemx 20
char A_ChNa[A_ShmShift+A_Nsemx][132]; /* Child process name (not used) */
int  A_ChID[A_ShmShift+A_Nsemx] = {0, 0, 0};        /* Child process ID */
int  A_ShmL[A_ShmShift+A_Nsemx]; /* Child Shmem segment length */
int  A_ShmID[A_ShmShift+A_Nsemx] = {0, 0};
void *A_ShmAdr[A_ShmShift+A_Nsemx];
            /* {sem_num, sem_op, sem_flag}; */
struct sembuf buf0 = {0, 0, ~SEM_UNDO&~IPC_NOWAIT};

/*-------------------------------
  Kill ASTRA, free memory
*/
void a_stop_(){
    printf("\nASTRA controlled stop\n");
    freeshm();
    int AID = getpid();
    printf("Kill process # %d\n", AID);
    kill(AID, SIGKILL);
    system("rm tmp/*.ppm >& /dev/null");
    exit(0);
}
/*---------------------------------------------------------------------*/
/* This function is called from Astra -> to_tra
   "to_tra" fills tra specific data in the dedicated memory segment
   Similarly, Astra -> ot_tra returns the data calculated
   by the process "tra" and stored in shared memory.
*/

void to_tra_(INT_* jrho_beg, INT_* jrho_end, INT_* N){
    struct shmid_ds Myshmid_ds;
    if (A_ShmNum < 0) return;
    if (shmctl(A_ShmID[*N+1], IPC_STAT, &Myshmid_ds) < 0){
        printf(">>> Process # %d: shmctl error >>>\n",*N+1);
        exit(1);
    }
    else if (Myshmid_ds.shm_segsz != A_ShmL[*N+1]){
        printf(">>> Process No.%d, \"%s\" >>>",*N,&A_ChNa[*N+1][0]);
        printf(" Size of shared memory mismatch\n");
        printf("    Allocated %d != %ld(found)\n",
            A_ShmL[*N+1], Myshmid_ds.shm_segsz);
        exit(1);
    }
    ql_io = (struct A_ql_io *)A_ShmAdr[*N+1];
    ql_io->jrho_beg = *jrho_beg;
    ql_io->jrho_end = *jrho_end;
    return;
}

/*----------------------------------------------------------------*/
void ot_tra_(INT_* jrho_beg, INT_* jrho_end, INT_* N, double* cpuse, double* mem){
// Reads the shared memory segment and stores the QL-code output to ASTRA fortran arrays 
    int j, jarr;
    if (A_ShmNum < 0) return;
    size_t offset = offsetof(struct A_ql_io, QLarrays);
    ql_io = (struct A_ql_io *)A_ShmAdr[*N+1];
    ql_io->QLarrays = (double *)((char *)ql_io + offset);
    *cpuse = ql_io->My.CPUse;
    for (j=*jrho_beg-1; j <= *jrho_end-1; j++){
        for (jarr=0; jarr<N_QL; jarr++){
            mem[j+1+jarr*A_NA1] = ql_io->QLarrays[j+jarr*A_NA1];
        }
    }
    return;
}

/*----------------------------------------------------------------*/
char* parse_nml(char * line_in){
    if (line_in == NULL) return NULL;

// Find the opening quote
    const char *start = strchr(line_in, '"');
    if (!start) start = strchr(line_in, '\'');
    if (!start) return NULL;

    start++;

// Find the end of the quoted word
    const char *end = start;
    while (*end && *end != '"' && *end != '\'' && !isspace((unsigned char)*end)) {
        end++;
    }

    size_t len = end - start;
    char *word_out = malloc(len + 1);  // +1 for null terminator
    if (!word_out) return NULL;

    strncpy(word_out, start, len);
    word_out[len] = '\0';

    return word_out;
}

/*-----------------------------------------------------
  Check existence of executable files listed in subs
  Reads tmp/astra.log and fills external variables AWD, MOD, DATA
  *Nsub - total number of files_names/strings in subs,
  *Lstr - length of an element of the character ARRAY "subs",
          maximum length of the subprocess_name,
  *subs - character ARRAY, described in a calling Fortran routine as
          character*(*Lstr) ARRAY(max_length)
   each element includes a name of external process to be called
*/

int checkexec_(INT_* Nsub, INT_ *Lstr, char *subs){
    char stri[132], name[132], path[132];
    char *line;
    size_t len = 0;
    ssize_t read;
    int j, i;
    if (A_Nsems != 0) return(A_Nsems); /* Do check only once */
/* Check existence of subprocess executable files */
/* Read tmp/astra.log and store run info */
    FILE *A_LOG;
    A_LOG = fopen(A_log_file, "r");
    if (!A_LOG){
        printf("Cannot open Astra log file: \"%s\"\n", A_log_file);
        exit(0);
    }

    while ((read = getline(&line, &len, A_LOG)) != -1) {
        if (strstr(line, "AWD"     ) != NULL) AWD    = parse_nml(line);
        if (strstr(line, "equ_file") != NULL) equmod = parse_nml(line);
        if (strstr(line, "exp_file") != NULL) DATA   = parse_nml(line);
    }

    fclose(A_LOG);
    free(line);

    for (j = 0; j < *Nsub; j++) {
        char *sub = &subs[*Lstr * j];

        if (strlen(sub) == 0) {
            fprintf(stderr, " >>> Xroutine call string (#%d) is empty.\n", j + 1);
            exit(EXIT_FAILURE);
        }

        strncpy(path, sub, sizeof(path) - 1);
        path[sizeof(path) - 1] = '\0';

// Expand ~ to $HOME if necessary
        if (path[0] == '~') {
            const char *home = getenv("HOME");
            if (!home) {
                fprintf(stderr, " >>> Xroutine call string \"%s\" error:\n", sub);
                fprintf(stderr, " >>> Symbol \"~\" is not allowed.\n");
                exit(EXIT_FAILURE);
            }

            snprintf(stri, sizeof(stri), "%s%s", home, path + 1);
            strncpy(path, stri, sizeof(path) - 1);
            path[sizeof(path) - 1] = '\0';
        }

// Separate name and path
        char *slash = strrchr(path, '/');
        if (slash) {
            strncpy(name, slash + 1, sizeof(name) - 1);
            name[sizeof(name) - 1] = '\0';
            *(slash + 1) = '\0'; // Truncate path after last slash
        } else {
            strncpy(name, path, sizeof(name) - 1);
            name[sizeof(name) - 1] = '\0';
            path[0] = '\0';
        }

        if (strlen(path) > 62) {
            fprintf(stderr, " >>> Xroutine call string \"%s\" error:\n", sub);
            fprintf(stderr, " >>> Absolute path is too long.\n");
            exit(EXIT_FAILURE);
        }

// Check if executable exists
        snprintf(stri, sizeof(stri), "test -x %s%s", path, name);
        if (system(stri) == 0) {
            A_Nsems++;
        } else {
            fprintf(stderr, "The executable file \"%s%s\" (#%d) does not exist\n", path, name, j + 1);
            exit(j);
        }
    }

    A_Nsems++;
    return(A_Nsems);
}

/*---------------------------------------------------
  "call initipc(NA1)" is placed in init.inc
  Get PID and key for the Astra main process
  Create and initialize a set of A_Nsems semaphores
  Assign NA1 (= *Ngrid) to A_NA1
  Allocate two shared memory segments for Astra datasets
*/
int initipc_(INT_* Ngrid, INT_ *n_ql){
    int l, var_size, arr_size, is=0, ds, j, *k;
    FILE *A_PDF;
    char hostname[132];
    size_t namlen;
    time_t hold_time;
    static union semun Mysemun;

    if (A_Nsems == 0) return(0); /* Remove this line if ESC is enabled */
    if (A_NA1 != 0) return(0); /* Initialize only once */

/* Collecting data */
    A_PID = getpid();
/* Define the absolute path name of Astra executable ASTRA_task */
    strcpy(ASTRA_task, AWD);
    strcat(ASTRA_task, "bin/");
    strcat(ASTRA_task, equmod);
    strcat(ASTRA_task, ".exe");
    my_key = ftok( ASTRA_task, (int)A_PID);    /* Get System V IPC key */
    printf("ASTRA_task %s\n", ASTRA_task);
    if (my_key == -1){
        printf("ipc_control: not able to create Key from ProcID\n");
        printf("Probably wrong ATASK name parsed from tmp/astra.nml\n");
        printf("ATASK name: %s\n", ASTRA_task);
        exit(1);
    }

    if (A_Nsems != 0){
/* Create a set of A_Nsems semaphores */
        A_SemID = semget(my_key, A_Nsems, 0660|IPC_CREAT);
/* Initialize all semaphores in the set A_SemID as {0, 0, ...} */
        Mysemun.val = 0;
        for(j=0; j < A_Nsems; j++){
            semctl(A_SemID, j, SETVAL, Mysemun);
        }
    }
/* Write file tmp/astra.ipc
   Form absolute path for A_ipc_file
*/
    strcpy(A_ipc_file, AWD);
    strcat(A_ipc_file, "/tmp/");
    strcat(A_ipc_file, DATA);
    strcat(A_ipc_file, equmod);
    strcat(A_ipc_file, ".ipc");
    A_PDF = fopen(A_ipc_file, "w");
    if (!A_PDF){
        printf("Cannot open Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(0);
    }
    fprintf(A_PDF, " Astra task:  \"%s\"\n", ASTRA_task);
    fprintf(A_PDF, " Astra files:  \"%s\",  \"%s\"\n", DATA, equmod);
    gethostname(hostname, (size_t)32);

    hold_time=time(NULL);
    fprintf(A_PDF, " Astra@%s started on:  %s", hostname, ctime(&hold_time));
/* The next line will be used if ESC is enabled */
    if (A_Nsems == 0){
        fclose(A_PDF);
        return(0);
    }
    fprintf(A_PDF, " Astra(main):  PID = %d,  SemID = %d\n",
    (int)A_PID, A_SemID);
    A_NA1 = *Ngrid;
    N_QL = *n_ql;

    var_size = sizeof(struct A_vars);
    AllocateShmem(var_size);
    fprintf(A_PDF, " Astra_inout: ShmID(A_vars):%12d%12d\n", A_ShmID[0], var_size);

#include "A_arrs.h"
    arr_size = sizeof(struct A_arrs);
    AllocateShmem(arr_size);
    fprintf(A_PDF, " Astra_inout: ShmID(A_arrs):%12d%12d\n", A_ShmID[1], arr_size);
    fprintf(A_PDF, "        PID      ShmID       ShmSize     Process\n");
    fclose(A_PDF);

    AVARS = (struct A_vars *)A_ShmAdr[0];
    AARRS = (struct A_arrs *)A_ShmAdr[1];
}

/*----------------------------------------------------------
  Get ShMemIDs for Astra datasets (const.inc) and (status.inc)
*/
void AllocateShmem(int l){
    if (A_ShmNum > A_ShmShift+A_Nsemx){
        printf(">>> ERROR >>> Too many shared memory segments requested\n");
        a_stop_();
    }
    A_ShmNum++;
/* Allocate a shared memory segment starting at A_ShmAdr[A_ShmNum] */
    A_ShmID[A_ShmNum] = shmget((key_t)(my_key+A_ShmNum), l, 0660|IPC_CREAT|IPC_EXCL);
/* Attach shared memory to the process */
    A_ShmAdr[A_ShmNum] = shmat(A_ShmID[A_ShmNum], NULL, 0);

    return;
}

/*---------------------------------------------------------------------
  Set (lock) the primary semaphore to -(Number_of_processes)
  Launch subordinate processes
  Those in turn do the following
    (1) get all PIDs, Keys, ShMems
    (2) create own ShMem segment
    (3) append file tmp/astra.ipc
    (4) increment primary semaphore
    (5) wait until the primary opens track
*/
int inikids_(INT_* Nsub, INT_ *Lstr, char *subs){
    if (A_ChID[A_ShmShift] != 0) return(0); /* Initialize only once */
    char stri[132], name[132], path[132];
    int j, i;

    if (A_Nsems <= *Nsub){
        printf(" >>> ERROR >>> Incomplete semaphore set\n");
       return(1);
    }

    for (j = 0; j < *Nsub; j++) {
        char *sub = &subs[*Lstr * j];

        if (strlen(sub) == 0) {
            printf("Error in input SBP string [%s]\n", sub);
            return j;
        }

        strncpy(path, sub, sizeof(path) - 1);
        path[sizeof(path) - 1] = '\0';

        const char *home = getenv("HOME");
        if (path[0] == '~' && home) {
            snprintf(stri, sizeof(stri), "%s%s", home, path + 1);
            strncpy(path, stri, sizeof(path) - 1);
            path[sizeof(path) - 1] = '\0';
        }

        char *slash = strrchr(path, '/');
        if (slash) {
            strncpy(name, slash + 1, sizeof(name) - 1);
            name[sizeof(name) - 1] = '\0';
            *(slash + 1) = '\0'; // Truncate path after last slash
        } else {
            strncpy(name, path, sizeof(name) - 1);
            name[sizeof(name) - 1] = '\0';
            path[0] = '\0';
        }

        int i;
        if (strlen(path) == 0) {
            snprintf(stri, sizeof(stri), "./bin/%s %s %d %d %d &",
                     name, ASTRA_task, A_PID, (int)my_key, j + 1);
            i = system(stri);
        } else {
            if (chdir(path) != 0) {
                perror("chdir failed");
                return j + 1;
            }

            snprintf(stri, sizeof(stri), "./%s %s %d %d %d &",
                     name, ASTRA_task, A_PID, (int)my_key, j + 1);
            i = system(stri);

            if (chdir(AWD) != 0) {
                perror("chdir back to AWD failed");
                return j + 1;
            }
        }

// Wait until child increments semaphore 0
        struct sembuf bufj = {0, -1, ~IPC_NOWAIT};
        if (semop(A_SemID, &bufj, 1) == -1) {
            perror("semop failed");
            return j + 1;
        }

        if (i == -1) return j + 1;
    }

/* All secondaries all launched and the file A_ipc_file is completed
   Now the data from A_ipc_file have to be retrieved by the main */
    if (read_aipc(Nsub, Lstr, subs)){
        printf("Error in input data interpretation (function read_aipc)\n");
        exit(j);
    }
    return(0);
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
    for (j=0; j < A_Nsems-1; j++) printf("%d, ", semarray[j]);
    printf("%d}\n", semarray[A_Nsems-1]);
    return;
}

/*--------------------- Check if ipc is activated --------------------*/
int ifipc_(){
    if (A_NA1 == 0) return(0);
    return(1);
}

/*--------------------- Unlock subprocess ----------------------------*/
int letsbp_(INT_* n){
    auto struct sembuf bufN = {*n, 1, IPC_NOWAIT};
    --buf0.sem_op;   /* Each call decrements Sem0 value by 1 */
/* Increment semval # sem_num=*n, Open subprocess */
    semop(A_SemID, &bufN, 1);
}

/*---------------------------------------------------------------------*/
/* Compare Sem0 value with buf0.sem_op ( == -A_Nsems) and
 wait until all subprocesses increment Sem0 by 1
        so that Sem0 reaches value A_Nsems                       ------*/
int wait4all_(){
    if (A_ShmNum < 0) return(0); /* Do check only after initialization   */
    static struct timespec timeout = {0, 100000000};   /* timeout = .1 sec */

    MinorLoop:

/*
  Go on if [Sem0value+buf0.sem_op==0], goto Minorloop after timeout
  Here: buf0 = {(ushort_t) buf0.sem_num = 0 = Const,
                -1 <= (short) buf0.sem_op <= Number_of_active_subprocesses,
        (short) buf0.sem_flg = SEM_UNDO = Const)}
  Note! Each call increments value of semadj by 1.
  Overflow occurs when the number of successive calls exceeds 32767
*/

    if (semtimedop(A_SemID, &buf0, 1, &timeout)){
        switch(errno){
        case EAGAIN:
            goto MinorLoop;
        case EIDRM:
            printf("The semaphore set was removed\n");
        case EINTR:
            printf("While blocked in this call, the process caught a signal\n");
        case ERANGE:
            printf("\n");
            WhatSem();
            printf("ERANGE error:  sem_op = %d,   SEMVMX = ?\n\n", buf0.sem_op);
        case EINVAL:
            printf("EINVAL:\t");
        case EFAULT:
            printf("EFAULT:\t");
        case E2BIG:
            printf("E2BIG:\t");
        case EACCES:
            printf("EACCES:\t");
        case EFBIG:
            printf("EFBIG:\t");
        case ENOMEM:
            printf("ENOMEM:\t");
        default:
            printf(">>> PRIMARY >>> Unrecognised semtimedop error\n");
        }
        printf(">>> PRIMARY >>> semtimedop error = %d\n", errno);
        a_stop_();
    }
    else{
        buf0.sem_op = 0;
        return(0);
    }

}

/*------------------------------------------------------------
  Read IPC file after launching all processes and
    (1) fill arrays ofChild_process_IDs, Child_shmem_lengths/IDs,
    (2) attach child process shmem segments to the main process memory
*/
int read_aipc(INT_* Nsub, INT_ *Lstr, char *subs){
    FILE *A_PDF;
    char stri[132], name[132];
    int j, i, k, ID, ShmID, kS;

    A_PDF = fopen(A_ipc_file, "r");
    if (!A_PDF){
        printf("Cannot open Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(0);
    }
    for (j=0; j < 5+A_ShmShift; j++){
        fgets(stri, 132, A_PDF); /* Skip lines */
    }
    i = A_ShmShift;
    while (EOF != fscanf(A_PDF, "%12d%12d%12d%s", &ID, &ShmID, &kS, stri) ){
        for (j = i-A_ShmShift; j < A_Nsems; j++){
            if ( strrchr(&subs[*Lstr*j], '/') != NULL ){
                k = strrchr(&subs[*Lstr*j], '/') - &subs[*Lstr*j];
                strcpy(name, &subs[*Lstr*j+k+1]);
            }
            else{
                strcpy(name, &subs[*Lstr*j]);
            }
            if (j == i-A_ShmShift){
                break;      /* SBP name matches the record */
            }
            else{
                goto Out;  /* Inconsistency in file A_ipc_file */
            }
        }
        if (A_ShmNum <= A_ShmShift+A_Nsemx){
            A_ShmNum++;
            A_ChID[A_ShmNum] = ID;
            A_ShmL[A_ShmNum] = kS;
            A_ShmID[A_ShmNum] = ShmID;    /* Here (A_ShmNum == Shmem_index) */
            A_ShmAdr[A_ShmNum] = shmat(ShmID, NULL, 0);
            strcpy(&A_ChNa[A_ShmNum][0], stri);
        }
        else{
            printf(">>> ERROR >>> Too many shmem segments\n");
            a_stop_();
        }
        i++;
    }
    fclose(A_PDF);
    if ( A_ShmNum+1-A_ShmShift == *Nsub ) return(0);

    Out:
    printf(" >>> File \"%s\" error >>> Missing SBP(%d) name: [%s]\n",
    A_ipc_file, i+1-A_ShmShift, name);
    printf("\n  File \"%s\" contents:\n", A_ipc_file);
    printf("\n  File processing:\n");
    printf("A_Nsems = %d,  A_ShmShift = %d,  A_ShmNum = %d,  Nsub = %d\n",
    A_Nsems, A_ShmShift, A_ShmNum, *Nsub);
    for (j=A_ShmShift-1; j <= *Nsub; j++){
        fprintf(stdout, "%12d%12d%12d%12d  %s\n",
        j-A_ShmShift+1, A_ChID[j], A_ShmL[j], A_ShmID[j],
        A_ChNa[j]);
    }
    return(1);
}

/*----------- Get ID of ShMem for ASTRA scalars -------------*/
/* First active only after "initipc", i.e. after "init.inc" */
int setvars_(double* DEVAR, INT_* NA1, INT_* NB1, INT_* NBOUND, INT_* N, INT_* N_QL){
    int *I, j;

    if (A_Nsems == 0){
        printf(" >>> setvars >>> Illegal call: Semaphores are not created\n");
        return(0);
    }
    if (A_ShmNum < 0) return(0);

    AVARS = (struct A_vars *)A_ShmAdr[0];
    AVARS->ab    = *(DEVAR);
    AVARS->abc   = *(DEVAR + 1);
    AVARS->aim1  = *(DEVAR + 2);
    AVARS->aim2  = *(DEVAR + 3);
    AVARS->aim3  = *(DEVAR + 4);
    AVARS->amj   = *(DEVAR + 5);
    AVARS->btor  = *(DEVAR + 7);
    AVARS->elong = *(DEVAR + 8);
    AVARS->encl  = *(DEVAR + 10);
    AVARS->enwm  = *(DEVAR + 11);
    AVARS->ipl   = *(DEVAR + 18);
    AVARS->nncl  = *(DEVAR + 20);
    AVARS->nnwm  = *(DEVAR + 21);
    AVARS->rtor  = *(DEVAR + 27);
    AVARS->shift = *(DEVAR + 28);
    AVARS->trian = *(DEVAR + 29);
    AVARS->updwn = *(DEVAR + 32);
    AVARS->zmj   = *(DEVAR + 36);
    AVARS->n_ql = *N_QL;
    AVARS->na1  = *NA1;
    AVARS->nb1  = *NB1;
    AVARS->nrd  = *N;
    AVARS->na1n = *(NBOUND);
    AVARS->na1e = *(NBOUND + 1);
    AVARS->na1i = *(NBOUND + 2);
    return(0);
}

/*------------- Get ID of ShMem for status.inc --------------*/
int setarrs_(double* plasma_profs, INT_ *n_RD){
    int jrho;
    int NRD = *n_RD;
    if (A_Nsems == 0) return(0);
    if (A_ShmNum < 0) return(0);
#include "A_arrs.h"
    AARRS = (struct A_arrs *)A_ShmAdr[1];
    for (jrho=0; jrho < A_NA1; jrho++){
        AARRS->amain[jrho] = *(plasma_profs + jrho);
        AARRS->ametr[jrho] = *(plasma_profs + jrho +    NRD);
        AARRS->cu[jrho]    = *(plasma_profs + jrho +  2*NRD);
        AARRS->elon[jrho]  = *(plasma_profs + jrho +  3*NRD);
        AARRS->er[jrho]    = *(plasma_profs + jrho +  4*NRD);
        AARRS->fp[jrho]    = *(plasma_profs + jrho +  5*NRD);
        AARRS->g11[jrho]   = *(plasma_profs + jrho +  6*NRD);
        AARRS->ipol[jrho]  = *(plasma_profs + jrho +  7*NRD);
        AARRS->mu[jrho]    = *(plasma_profs + jrho +  8*NRD);
        AARRS->nalf[jrho]  = *(plasma_profs + jrho +  9*NRD);
        AARRS->ndeut[jrho] = *(plasma_profs + jrho + 10*NRD);
        AARRS->ne[jrho]    = *(plasma_profs + jrho + 11*NRD);
        AARRS->nhe3[jrho]  = *(plasma_profs + jrho + 12*NRD);
        AARRS->nhydr[jrho] = *(plasma_profs + jrho + 13*NRD);
        AARRS->ni[jrho]    = *(plasma_profs + jrho + 14*NRD);
        AARRS->nibm[jrho]  = *(plasma_profs + jrho + 15*NRD);
        AARRS->niz1[jrho]  = *(plasma_profs + jrho + 16*NRD);
        AARRS->niz2[jrho]  = *(plasma_profs + jrho + 17*NRD);
        AARRS->niz3[jrho]  = *(plasma_profs + jrho + 18*NRD);
        AARRS->ntrit[jrho] = *(plasma_profs + jrho + 19*NRD);
        AARRS->pblon[jrho] = *(plasma_profs + jrho + 20*NRD);
        AARRS->pbper[jrho] = *(plasma_profs + jrho + 21*NRD);
        AARRS->pfast[jrho] = *(plasma_profs + jrho + 22*NRD);
        AARRS->rho[jrho]   = *(plasma_profs + jrho + 23*NRD);
        AARRS->shear[jrho] = *(plasma_profs + jrho + 24*NRD);
        AARRS->shif[jrho]  = *(plasma_profs + jrho + 25*NRD);
        AARRS->te[jrho]    = *(plasma_profs + jrho + 26*NRD);
        AARRS->ti[jrho]    = *(plasma_profs + jrho + 27*NRD);
        AARRS->tria[jrho]  = *(plasma_profs + jrho + 28*NRD);
        AARRS->upl[jrho]   = *(plasma_profs + jrho + 29*NRD);
        AARRS->vpol[jrho]  = *(plasma_profs + jrho + 30*NRD);
        AARRS->vrs[jrho]   = *(plasma_profs + jrho + 31*NRD);
        AARRS->vtor[jrho]  = *(plasma_profs + jrho + 32*NRD);
        AARRS->zef[jrho]   = *(plasma_profs + jrho + 33*NRD);
        AARRS->zim1[jrho]  = *(plasma_profs + jrho + 34*NRD);
        AARRS->zim2[jrho]  = *(plasma_profs + jrho + 35*NRD);
        AARRS->zim3[jrho]  = *(plasma_profs + jrho + 36*NRD);
        AARRS->zmain[jrho] = *(plasma_profs + jrho + 37*NRD);
        AARRS->zimpt[jrho] = *(plasma_profs + jrho + 38*NRD);
        AARRS->nimpt[jrho] = *(plasma_profs + jrho + 39*NRD);
        AARRS->aimpt[jrho] = *(plasma_profs + jrho + 40*NRD);
    }
    return(0);
}

/*---------------------------------------------------------------------*/
void freeshm(){
    struct shmid_ds Myshmid_ds;
    int j;
    printf(">>> Freeing shared memory segments\n");
    if (A_Nsems != 0) semctl(A_SemID, 0, IPC_RMID);
    if (A_ShmNum < 0) return;
    for (j=0; j <= A_ShmNum; j++){
/* Detach and remove all shared memory segments */
        if (shmdt(A_ShmAdr[j]) < 0){
            if (errno == EINVAL){
                printf(">>>  Invalid ShmAdr[j] value: %d\n", j);
            }
        }
        else if (shmctl(A_ShmID[j], IPC_RMID, &Myshmid_ds) < 0){
            if (errno == EPERM){
                printf(">>> User cannot remove shared memory segment #%d\n", j);
                break;
            }
            else{
                printf(">>> Process # %d >>>  unknown shmctl error\n", j);
            }
        }
    }
    return;
}

/*---------------------------------------------------------------------*/
void write_aipc(const struct A_proc_info Aproc, char* AWD, int qlSize)
{
    FILE *A_PDF, *A_LOG;
    char A_logf[132];
    char *line;
    size_t len=0;
    ssize_t read;

/* Check existence of subprocess executable files */
/* Read tmp/astra.nml and store run info */

    strcpy(A_logf, AWD);
    strcat(A_logf, "/tmp/astra.nml");
    A_LOG = fopen(A_logf, "r");
    if (!A_LOG){
        printf("Cannot open Astra log file: \"%s\"\n", A_log_file);
        exit(0);
    }

    while ((read = getline(&line, &len, A_LOG)) != -1) {
        if (strstr(line, "equ_file") != NULL) equmod = parse_nml(line);
        if (strstr(line, "exp_file") != NULL) DATA = parse_nml(line);
    }

    fclose(A_LOG);
    free(line);

    strcpy(A_ipc_file, AWD);
    strcat(A_ipc_file, "tmp/");
    strcat(A_ipc_file, DATA);
    strcat(A_ipc_file, equmod);
    strcat(A_ipc_file, ".ipc");
    A_PDF = fopen(A_ipc_file, "a");
    if (!A_PDF){
        printf("Cannot open existing Astra IPC file: \"%s\"\n", A_ipc_file);
        exit(0);
    }
    fprintf(A_PDF, "%12d%12d%12d   %s\n", getpid(), Aproc.ShMid, qlSize, Aproc.Path);
    fclose(A_PDF);
    return;
}

/*---------------------------------------------------------------------*/
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
    }  // Overflow range of clock_t
    if (time0 == 0){     // Set time0 at start
        secs_per_tick = 1./sysconf(_SC_CLK_TCK);
        time0 = times(&buf);
        prev_time = buf.tms_utime + buf.tms_stime;
        return 0.;
    }
    run_time = times(&buf) - time0;  // Set time difference
    cpu_time = buf.tms_utime + buf.tms_stime;
    *secs += (cpu_time - prev_time)*secs_per_tick;
    prev_time = cpu_time; 
    runsec = run_time*secs_per_tick;
    return runsec;
}
