#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <signal.h>
#include <ctype.h>
#include <stddef.h>

#ifndef INT8
#define INT_ int
#else
#define INT_ long
#endif
#define IPCACTIVE 1

#include "A_vars.h"
#include "A_proc.h"

void AstraEvent();
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

union semun{
    int val;                    /* value for SETVAL */
    struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
    unsigned short int *array;  /* array for GETALL, SETALL */
    struct seminfo *__buf;      /* buffer for IPC_INFO */
};

pid_t A_PID = 0;
INT_  A_NB1 = 0;
#define NC1 A_NB1
int A_SemID = 0;
int A_Nsems = 0;       /* the number of semaphores */
int A_ShmNum = -1;
#define A_ShmShift 2
#define A_Nsemx    20
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
#ifdef IPCACTIVE
    freeshm();
#endif
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
    static struct A_ql_in
    {
        struct A_proc_info My; /* General IO information */
        int Size;         /* Control: Size of the Shmem */
        int jrho_beg;
        int jrho_end;
    } *QL_IN;
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
    QL_IN = (struct A_ql_in *)A_ShmAdr[*N+1];
    (*QL_IN).jrho_beg = *jrho_beg;
    (*QL_IN).jrho_end = *jrho_end;
    return;
}

/*----------------------------------------------------------------*/
void ot_tra_(INT_* jrho_beg, INT_* jrho_end, INT_* N, double* cpuse, double* YY){

#include "A_ql_io2.h"
    int j, jarr, n_nrd;
    if (A_ShmNum < 0) return;
    AVARS = (struct A_vars *)A_ShmAdr[0];
    n_nrd = AVARS->nrd;
    size_t offset = offsetof(struct A_ql_io2, QLarrays);
    ql_io2 = (struct A_ql_io2 *)A_ShmAdr[*N+1];
    ql_io2->QLarrays = (double *)((char *)ql_io2 + offset);
    *cpuse = ql_io2->My.CPUse;
    for (j=*jrho_beg-1; j <= *jrho_end-1; j++){
        for (jarr=0; jarr<15; jarr++){
            YY[j+1+jarr*n_nrd] = ql_io2->QLarrays[j+jarr*NC1];
        }
    }
    if (*jrho_beg == 1){
        for (j=0; j <= 15*n_nrd; j += n_nrd) YY[j] = 0.;
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

    start++;  // Move past opening quote

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

    for(j=0; j < *Nsub; j++){
        if (strlen(&subs[*Lstr*j]) == 0) goto Error1;
        strcpy(path, &subs[*Lstr*j]);

        if (strchr(path, '~') != NULL){
            if ( getenv("HOME") == NULL ) goto Error2;
            strcpy(stri, getenv("HOME"));
            strcat(stri, &path[1]);
            strcpy(path, stri);
            i = strrchr(path, '/') - &path[0];
            strcpy(name, &path[i+1]);
            if (strlen(path) > 62) goto Error3;
            path[i+1] = '\0';
        }
        else if ( strrchr(path, '/') != NULL ){
            i = strrchr(path, '/') - &path[0];
            path[i+1] = '\0';
            strcpy(name, &subs[*Lstr*j+i+1]);
        }
        else{
            strcpy(name, path);
            path[0] = '\0';
        }

        strcpy(stri, "test -x "); // Check for xpr/tglfi executable
        strcat(stri, path);
        strcat(stri, name);

        if (system(stri) == 0){
            A_Nsems++;
        }
        else{
            printf("The executable file \"%s%s\" (#%d) does not exist\n", 
                path, name, j+1);
            exit(j);
        }
    }
    A_Nsems++;
    return(A_Nsems);

    Error1:
    printf(" >>> Xroutine call string (#%d) \"%s\" error\n", j+1, &subs[*Lstr*j]);
    exit(0);
    Error2:
    printf(" >>> Xroutine call string \"%s\" error:\n", &subs[*Lstr*j]);
    printf(" >>> Symbol \"~\" is not allowed.\n");
    exit(0);
    Error3:
    printf(" >>> Xroutine call string \"%s\" error:\n", &subs[*Lstr*j]);
    printf(" >>> Absolute path is too long.\n");
    exit(0);
}

/*---------------------------------------------------
  "call initipc(NB1)" is placed in init.inc
  Get PID and key for the Astra main process
  Create and initialize a set of A_Nsems semaphores
  Assign NB1 (= *Ngrid) to A_NB1 (alias NC1)
  Allocate two shared memory segments for Astra datasets 
*/
int initipc_(INT_* Ngrid){
    int l, var_size, arr_size, is=0, ds, j, *k;
    FILE *A_PDF;
    char hostname[132];
    size_t namlen;
    time_t hold_time;
    static union semun Mysemun;

    if (A_Nsems == 0) return(0); /* Remove this line if ESC is enabled */
    if (A_NB1 != 0) return(0); /* Initialize only once */

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
    A_NB1 = *Ngrid;

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
    AVARS->c_size[0] = var_size;
    AVARS->c_size[1] = sizeof(int);
    AVARS->c_size[2] = sizeof(double);
    AVARS->c_size[3] = sizeof(key_t);
    AVARS->c_size[4] = sizeof(pid_t);
    AVARS->CheckWord = 314159265;

    AARRS = (struct A_arrs *)A_ShmAdr[1];
    AARRS->NC1  = NC1;
    AARRS->Size = arr_size;
    AARRS->CheckWord = 314159265;
}

/*----------------------------------------------------------
  Get ShMemIDs for Astra datasets (const.inc) and (status.inc)
*/
void AllocateShmem (int l){
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

    for(j=0; j < *Nsub; j++){
        if (strlen(&subs[*Lstr*j]) == 0) goto Error1;
        strcpy(path, &subs[*Lstr*j]);
        if (strchr(path, '~') != NULL ){
            strcpy(stri, getenv("HOME"));
            strcat(stri, &path[1]);
            strcpy(path, stri);
            i = strrchr(path, '/') - &path[0];
            strcpy(name, &path[i+1]);
            path[i+1] = '\0';
        }
        else if (strrchr(path, '/') != NULL){
            i = strrchr(path, '/') - &path[0];
            path[i+1] = '\0';
            strcpy(name, &subs[*Lstr*j+i+1]);
        }
        else{
            strcpy(name, path);
            path[0] = '\0';
        }
        i = strlen(path);

        if ( i == 0 ){
            sprintf(stri, "%s%s %s %d %d %d &", 
            "./bin/", name, ASTRA_task, A_PID, (int)my_key, j+1);
            i = system(stri);
        }
        else{
// Note! chdir() does not recognize ~ as home directory
            chdir(path); 
            sprintf(stri, "./%s %s %d %d %d &", 
            name, ASTRA_task, A_PID, (int)my_key, j+1);
            i = system(stri);
            chdir(AWD);
        }

/* Wait until proc No.(j+1) opens the PRIMARY semaphore (sem_num=0)
                            incrementing its initial value (0) by 1 */
        struct sembuf bufj = {0, -1, ~IPC_NOWAIT};
        semop(A_SemID, &bufj, 1);

        if (i == -1) return(j+1);
    }
/* All secondaries all launched and the file A_ipc_file is completed
   Now the data from A_ipc_file have to be retrieved by the main */
    if (read_aipc(Nsub, Lstr, subs)) goto Error2;
    return(0);

    Error1:
    printf("Error in input SBP string [%s]\n", &subs[*Lstr*j]);
    return(j);
    Error2:
    printf("Error in input data interpretation (function read_aipc)\n");
    exit(j);
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
    if (A_NB1 == 0) return(0);
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

    MinorLoop:{
/* Here the primary process can do limited actions e.g. analyze keys */
        (void) AstraEvent();
    }

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
                goto Out1;  /* Inconsistency in file A_ipc_file */
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

    Out1:
    printf(" >>> File \"%s\" error >>> Missing SBP(%d) name: [%s]\n", 
    A_ipc_file, i+1-A_ShmShift, name);
    goto Out;
    Out2:
    printf(" >>> File \"%s\" error >>> Incomplete record.\n", A_ipc_file);
    goto Out;
    Out:
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
int setvars_(double* DEVAR, INT_* NA1, INT_* NB1, INT_* NBOUND, INT_* N){
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
    AVARS->na1 = *NA1;
    AVARS->nb1 = *NB1;
    AVARS->nrd = *N;
    AVARS->na1n = *(NBOUND);
    AVARS->na1e = *(NBOUND + 1);
    AVARS->na1i = *(NBOUND + 2);
    return(0);
}

/*------------- Get ID of ShMem for status.inc --------------*/
int setarrs_(double* plasma_profs, INT_* NRD){
    int jrho;

    if (A_Nsems == 0) return(0);
    if (A_ShmNum < 0) return(0);
#include "A_arrs.h"
    AARRS = (struct A_arrs *)A_ShmAdr[1];
    for (jrho=0; jrho < NC1; jrho++){
        AARRS->amain[jrho] = *(plasma_profs + jrho); 
        AARRS->ametr[jrho] = *(plasma_profs + jrho + *NRD); 
        AARRS->cu[jrho]    = *(plasma_profs + jrho + 2*(*NRD)); 
        AARRS->elon[jrho]  = *(plasma_profs + jrho + 3*(*NRD)); 
        AARRS->er[jrho]    = *(plasma_profs + jrho + 4*(*NRD)); 
        AARRS->fp[jrho]    = *(plasma_profs + jrho + 5*(*NRD)); 
        AARRS->g11[jrho]   = *(plasma_profs + jrho + 6*(*NRD)); 
        AARRS->ipol[jrho]  = *(plasma_profs + jrho + 7*(*NRD)); 
        AARRS->mu[jrho]    = *(plasma_profs + jrho + 8*(*NRD)); 
        AARRS->nalf[jrho]  = *(plasma_profs + jrho + 9*(*NRD)); 
        AARRS->ndeut[jrho] = *(plasma_profs + jrho + 10*(*NRD)); 
        AARRS->ne[jrho]    = *(plasma_profs + jrho + 11*(*NRD)); 
        AARRS->nhe3[jrho]  = *(plasma_profs + jrho + 12*(*NRD)); 
        AARRS->nhydr[jrho] = *(plasma_profs + jrho + 13*(*NRD)); 
        AARRS->ni[jrho]    = *(plasma_profs + jrho + 14*(*NRD)); 
        AARRS->nibm[jrho]  = *(plasma_profs + jrho + 15*(*NRD)); 
        AARRS->niz1[jrho]  = *(plasma_profs + jrho + 16*(*NRD)); 
        AARRS->niz2[jrho]  = *(plasma_profs + jrho + 17*(*NRD)); 
        AARRS->niz3[jrho]  = *(plasma_profs + jrho + 18*(*NRD)); 
        AARRS->ntrit[jrho] = *(plasma_profs + jrho + 19*(*NRD)); 
        AARRS->pblon[jrho] = *(plasma_profs + jrho + 20*(*NRD)); 
        AARRS->pbper[jrho] = *(plasma_profs + jrho + 21*(*NRD)); 
        AARRS->pfast[jrho] = *(plasma_profs + jrho + 22*(*NRD)); 
        AARRS->rho[jrho]   = *(plasma_profs + jrho + 23*(*NRD)); 
        AARRS->shear[jrho] = *(plasma_profs + jrho + 24*(*NRD)); 
        AARRS->shif[jrho]  = *(plasma_profs + jrho + 25*(*NRD)); 
        AARRS->te[jrho]    = *(plasma_profs + jrho + 26*(*NRD)); 
        AARRS->ti[jrho]    = *(plasma_profs + jrho + 27*(*NRD)); 
        AARRS->tria[jrho]  = *(plasma_profs + jrho + 28*(*NRD)); 
        AARRS->upl[jrho]   = *(plasma_profs + jrho + 29*(*NRD)); 
        AARRS->vpol[jrho]  = *(plasma_profs + jrho + 30*(*NRD)); 
        AARRS->vrs[jrho]   = *(plasma_profs + jrho + 31*(*NRD)); 
        AARRS->vtor[jrho]  = *(plasma_profs + jrho + 32*(*NRD)); 
        AARRS->zef[jrho]   = *(plasma_profs + jrho + 33*(*NRD)); 
        AARRS->zim1[jrho]  = *(plasma_profs + jrho + 34*(*NRD)); 
        AARRS->zim2[jrho]  = *(plasma_profs + jrho + 35*(*NRD)); 
        AARRS->zim3[jrho]  = *(plasma_profs + jrho + 36*(*NRD)); 
        AARRS->zmain[jrho] = *(plasma_profs + jrho + 37*(*NRD)); 
        AARRS->zimpt[jrho] = *(plasma_profs + jrho + 38*(*NRD)); 
        AARRS->nimpt[jrho] = *(plasma_profs + jrho + 39*(*NRD)); 
        AARRS->aimpt[jrho] = *(plasma_profs + jrho + 40*(*NRD)); 
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

/*------- Call as:  write_aipc_(Aproc.OrdNr, Aproc.ShMid, lS), -------*/

void write_aipc (const struct A_proc_info Aproc, char* AWD, int* lS)
{
    FILE *A_PDF;
    char A_IPC[132]; 
    char A_logf[132];

    char *line;
    size_t len = 0;
    ssize_t read;
    int j, i;

/* Check existence of subprocess executable files */
/* Read tmp/astra.log and store run info */
    FILE *A_LOG;

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

    strcpy(A_IPC, AWD);
    strcat(A_IPC, "tmp/");
    strcat(A_IPC, DATA);
    strcat(A_IPC, equmod);
    strcat(A_IPC, ".ipc");
    A_PDF = fopen(A_IPC, "a");
    if (!A_PDF){
        printf("Cannot open existing Astra IPC file: \"%s\"\n", A_IPC);
        exit(0);
    }
    fprintf(A_PDF, "%12d%12d%12d   %s\n", getpid(), Aproc.ShMid, *lS, Aproc.Path);
    fclose(A_PDF);
    return;
}
/*-----------------------------------------------------
  Communication report
  Call as:  SP_stamp(Mama, Aproc, SemID);
*/
void SP_stamp(const struct A_proc_info Mama, 
              const struct A_proc_info Aproc, int SemID){

//    printf("#%d \"%s\", \tAproc_PID %d,  SemID %d\n  A_PID %d,  A_Key %d\n", 
//         Aproc.OrdNr, Aproc.Path, Aproc.Pid, SemID, Mama.Pid, Mama.Key); 
    return;
}
