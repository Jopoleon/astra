#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <errno.h>

typedef int32_t INT_;

static struct A_proc_info
{
  char   Path[96];      /* Path of a process */
  key_t  Key;           /* Key of a process */
  pid_t  Pid;           /* pid of a process */
  int    OrdNr;         /* Ordinal number of a process */
  int    ShMid;         /* ID of ShMem for a process */
  int    Size;          /* Size of A_proc_info */
  int    CheckWord;     /* Control number */
} A_proc_info;

static struct A_ql_io
{
   struct A_proc_info My; /* General IO information */
   int Size;         /* Control: Size of the Shmem */
   int jrho_beg;           /* Input:   */
   int jrho_end;           /* Input:   */
   double* QLarrays;
} *ql_io;

static struct A_vars
{
  double ab;
  double abc;
  double aim1;
  double aim2;
  double aim3;
  double amj;
  double btor;
  double elong;
  double encl;
  double enwm;
  double ipl;
  double nncl;
  double nnwm;
  double rtor;
  double shift;
  double trian;
  double updwn;
  double zmj;
  double time;
  double hro;
  double roc;
  double tau;
  int n_ql;
  int na1;
  int nab;
  int nb1;
  int nrd;
  int na1n;
  int na1e;
  int na1i; 
} *AVARS;

union semun{
    int val;                    /* value for SETVAL */
    struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
    unsigned short int *array;  /* array for GETALL, SETALL */
    struct seminfo *__buf;      /* buffer for IPC_INFO */
};
