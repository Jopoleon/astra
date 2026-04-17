#!/usr/bin/env python3

import os, sys, logging, argparse, traceback
import numpy as np
import greenMatrices

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')

logger = logging.getLogger('as_exe')

if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

f_log = '%s/tmp/astra.nml' %awd


def parse_alog():

    alog_d = {}
    if not os.path.isfile(f_log):
        return None
    for lin in open(f_log):
       keyval = lin.split('=')
       if len(keyval) > 1:
           key, val = keyval
           if val.strip() != '':
               alog_d[key.strip()] = val.strip().strip("'").strip('"')
    return alog_d


if __name__ == '__main__':

    sys.path.append(awd)
    import postProcess.json2cdf as j2nc

    greenMatrices.main()
    alog_d = parse_alog()
    parser = argparse.ArgumentParser(description='Write settings & run ASTRA')
    if alog_d is None:
        parser.add_argument('-m', '--equ', help='Model file', required=True)
        parser.add_argument('-v', '--exp', help='Exp file', required=True)
        parser.add_argument('-s', '--tbeg', type=float, help='Initial time', required=False, default=0.1)
        parser.add_argument('-e', '--tend', type=float, help='End time'    , required=False, default=10.)
        parser.add_argument('-dev', '--DeviceName', help='Machine short name', required=False, default='aug')
        parser.add_argument('-resize', '--resize', type=float, help='Resize factor for ASTRA frame', required=False, default=1.)
    else:
        parser.add_argument('-m', '--equ', help='Model file', required=False, default=alog_d['equ_file'])
        parser.add_argument('-v', '--exp', help='Exp file'  , required=False, default=alog_d['exp_file'])
        parser.add_argument('-s', '--tbeg', type=float, help='Initial time'  , required=False, default=alog_d['tbeg_nml'])
        parser.add_argument('-e', '--tend', type=float, help='End time'      , required=False, default=alog_d['tend_nml'])
        parser.add_argument('-dev', '--DeviceName', help='Machine short name', required=False, default=alog_d['machine'])
        parser.add_argument('-resize', '--resize', type=float, help='Resize factor for ASTRA frame', required=False, default=alog_d['resize'])
    parser.add_argument('-batch', action='store_true', help='Run batch job'  , required=False)
    parser.add_argument('-tpause', '--tpause', type=float, help='Pause @time', required=False, default=1.e4)
    parser.add_argument('-debug', action='store_true', help='Debug', required=False)
    parser.add_argument('-re', '--restart', type=int, help='Restart', required=False, default=0)
    parser.add_argument('-fs', action='store_true', help='Flight simulator', required=False)
    parser.add_argument('-W', '--waitslurm', action='store_true', help='Hold on SLURM job', required=False, default=False)
    parser.add_argument('-c', '--compiler', help='Compiler (gui, ifx)', required=False)

    args = parser.parse_args()

    dev_name = args.DeviceName
    resize = args.resize

    if args.batch:
        rtype = 'BGD'
    else:
        rtype = 'RUN'

    if args.compiler:
        os.environ['ASTRA_COMPILER'] = args.compiler

    if resize is None:
        if 'resize' in alog_d.keys():
            resize = float(alog_d['resize'])
        else:
            resize = 1
    if dev_name.strip() == '':
        dev_name = 'aug'

    alog  = '&astra_log\n\n'
    alog += 'exp_file  = "%s"\n'   %args.exp
    alog += 'equ_file  = "%s"\n'   %args.equ
    alog += 'tbeg_nml  = %8.4f\n'  %args.tbeg
    alog += 'tend_nml  = %8.4f\n'  %args.tend
    alog += 'TASK      = "%s"\n'   %rtype
    alog += 'machine   = "%s"\n'   %dev_name
    alog += 'debug     = %d\n'     %int(args.debug)
    alog += 'flightsim = .%s.\n'   %args.fs
    alog += 'resize    = %8.4f\n'  %resize
    alog += 'restart   = %d\n'     %args.restart
    alog += 'tpause_nml = %8.4f\n' %args.tpause
    alog += '\n/\n'

    logger.info('Writing %s' %f_log)
    with open(f_log, 'w') as f:
        f.write(alog)

    logger.debug(alog)

    expequ = args.exp + args.equ

    cmd = '%s/exe/Build %s %s %s %s' %(awd, args.equ, args.exp, rtype, args.fs)
    if args.batch and args.waitslurm:
        cmd += ' -W'
    logger.info(cmd)
    os.system(cmd)
# Remove pending IPC processes
    for ipcId in range(3):
        f_ipc = '%s/tmp/%s-%d.ipc' %(awd, expequ, ipcId)
        if os.path.isfile(f_ipc):
            with open(f_ipc) as f:
                lines = f.readlines()
            ipcProc1 = int(lines[6].split(':')[1].split()[0])
            ipcProc2 = int(lines[7].split(':')[1].split()[0])
            ipcProc3 = int(lines[8].split(':')[1].split()[0])
            jproc, PIDs, ipcProcs = np.loadtxt(f_ipc, skiprows=10, unpack=True, dtype=np.int32)
            for pid in PIDs:
                cmd = 'kill -9 %s' %pid
                os.system(cmd)
            for proc in np.append(ipcProcs, [ipcProc1, ipcProc2, ipcProc3]):
                cmd = 'ipcrm -m %d 2>/dev/null' %proc
                os.system(cmd)
    try:
        j2nc.json_concat(expequ)
    except:
#        logger.debug(traceback.format_exc())
        pass
