#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#ifndef INT8
#define INT_ int
#else
#define INT_ long
#endif

int to_tra_(INT_*, INT_*, INT_*);
int ot_tra_(INT_*, INT_*, INT_*, double*, double*);

extern char A_ChNa[][16]; /* Child process ID */
extern INT_ A_NB1;
extern int  A_ShmNum, A_ShmL[], A_ShmID[];
extern void *A_ShmAdr[];
#define NC1  A_NB1
#include "A_vars.h"
/*---------------------------------------------------------------------*/
/* This function is called from Astra -> totra -> to_tra
   "to_tra" fills tra specific data in the dedicated memory segment   
   Similarly, Astra -> ottra -> ot_tra returns the data calculated
   by the process "tra" and stored in shared memory.
*/

int to_tra_(INT_* IS, INT_* IE, INT_* N){

    int i, j, k;
    struct shmid_ds Myshmid_ds;
#include "A_proc.h"
#include "A_ql_IO.h"
    if (A_ShmNum < 0) return(0);
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

    AVARS = (struct A_vars *)A_ShmAdr[0];
    k = (*AVARS).nrd;
    IOQL = (struct A_ql_IO *)A_ShmAdr[*N+1];
    (*IOQL).is = *IS;
    (*IOQL).ie = *IE;
}

/*----------------------------------------------------------------*/
int ot_tra_(INT_* IS, INT_* IE, INT_* N, double* cpuse, double* YY){
    int j, i, n_nrd;
#include "A_proc.h"
#include "A_ql_IO.h"
    if (A_ShmNum < 0) return(0);
    AVARS = (struct A_vars *)A_ShmAdr[0];
    IOQL = (struct A_ql_IO *)A_ShmAdr[*N+1];
    *cpuse = (*IOQL).My.CPUse;
    n_nrd = AVARS->nrd;
    for (j=*IS-1; j <= *IE-1; j++){
        i = 1;
        YY[j+i] = (*IOQL).chi[j];  i += n_nrd; // work(j+1,1)
        YY[j+i] = (*IOQL).che[j];  i += n_nrd; // work(j+1,2)
        YY[j+i] = (*IOQL).dif[j];  i += n_nrd; // work(j+1,3)
        YY[j+i] = (*IOQL).vin[j];  i += n_nrd; // work(j+1,4)
        YY[j+i] = (*IOQL).dph[j];  i += n_nrd; // work(j+1,5)
        YY[j+i] = (*IOQL).dpl[j];  i += n_nrd; // work(j+1,6)
        YY[j+i] = (*IOQL).dpr[j];  i += n_nrd; // work(j+1,7)
        YY[j+i] = (*IOQL).xtb[j];  i += n_nrd; // work(j+1,8)
        YY[j+i] = (*IOQL).egm[j];  i += n_nrd; // work(j+1,9)
        YY[j+i] = (*IOQL).gam[j];  i += n_nrd; // work(j+1,10)
        YY[j+i] = (*IOQL).gm1[j];  i += n_nrd; // work(j+1,11)
        YY[j+i] = (*IOQL).gm2[j];  i += n_nrd; // work(j+1,12)
        YY[j+i] = (*IOQL).om1[j];  i += n_nrd; // work(j+1,13)
        YY[j+i] = (*IOQL).om2[j];  i += n_nrd; // work(j+1,14)
        YY[j+i] = (*IOQL).fr1[j];  i += n_nrd; // work(j+1,14)
    }
    if (*IS == 1){
       for (j=0; j <= 15*n_nrd; j += n_nrd) YY[j] = 0.;
    }
}
