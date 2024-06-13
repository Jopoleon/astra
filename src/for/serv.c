#include <unistd.h>
#include <stdio.h>
#include <sys/times.h>

/**********************************************************************/
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
    }  /* Overflow range of clock_t*/
    if (time0 ==  0){     /* Set time0 at start */
        secs_per_tick = 1./sysconf(_SC_CLK_TCK);
        time0 = times(&buf);
        prev_time = buf.tms_utime+buf.tms_stime;
        return 0.;
    }
    run_time = times(&buf) - time0;  /* Set time difference */
    cpu_time = buf.tms_utime + buf.tms_stime;
    *secs += (cpu_time - prev_time)*secs_per_tick;
    prev_time = cpu_time; 
    runsec = run_time*secs_per_tick;
    return runsec;
}
