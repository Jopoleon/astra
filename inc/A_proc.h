static struct A_proc_info
{
  char	 Path[96];      /* Path of a process */
  key_t	 Key;           /* Key of a process */
  pid_t	 Pid;           /* pid of a process */
  int	 OrdNr;         /* Ordinal number of a process */
  int	 ShMid;         /* ID of ShMem for a process */
  double CPUse;         /* CPU time uzage of a process */
  int	 Size;          /* Size of A_proc_info */
  int    CheckWord;     /* Control number */
} A_proc_info;
