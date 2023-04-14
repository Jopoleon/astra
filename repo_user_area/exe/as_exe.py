import os, logging, argparse, shutil, traceback
from nc_concat import nc_concat

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')

logger = logging.getLogger('as_exe')

if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

awd = os.getenv('AWD')
if awd is None:
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

    alog_d = parse_alog()
    parser = argparse.ArgumentParser(description='astra.nml writer')
    parser.add_argument('-m', '--equ', help='Model file', required=False, default=alog_d['equ_file'])
    parser.add_argument('-v', '--exp', help='Exp file', required=False, default=alog_d['exp_file'])
    parser.add_argument('-s', '--tbeg', type=float, help='Initial time', required=False, default=alog_d['tbeg_nml'])
    parser.add_argument('-e', '--tend', type=float, help='End time', required=False, default=alog_d['tend_nml'])
    parser.add_argument('-dev', '--DeviceName', help='Machine short name', required=False, default=alog_d['machine'])
    parser.add_argument('-batch', action='store_true', help='Run batch job', required=False)
    parser.add_argument('-tpause', '--tpause', type=float, help='Pause @time', required=False)
    parser.add_argument('-resize', '--resize', type=float, help='Resize factor for ASTRA frame', required=False)
    parser.add_argument('-debug', action='store_true', help='Debug', required=False)
    parser.add_argument('-fs', action='store_true', help='Flight simulator', required=False)

    args = parser.parse_args()

    dev_name = args.DeviceName
    resize = args.resize

    if args.batch:
        rtype = 'BGD'
    else:
        rtype = 'RUN'

    if resize is None:
        if 'resize' in alog_d.keys():
            resize = float(alog_d['resize'])
        else:
            resize = 1
    if dev_name.strip() == '':
        dev_name = 'aug'

    if '/afs/ipp' in awd:
        aext = '/afs/ipp/home/a/astra/ASTRA_LIBRARIES_EXT'
    elif '/toks/work' in awd:
        aext = '/afs/ipp/home/a/astra/ASTRA_LIBRARIES_EXT'
    elif '/home/ITER' in awd:
        aext = '/home/ITER/tarding/ASTRA_LIBRARIES_EXT'
    elif 'compass' in awd:
        aext = '/compass/home/tardini/astra_ext'
    elif 'eufus' in awd:
        aext = '/afs/eufus.eu/g2itmdev/user/g2gtardi/ASTRA_LIBRARIES_EXT'
    else:
        aext = '/fusion/projects/codes/astra/ASTRA_LIBRARIES_EXT'
    logger.debug(aext)
    alog  = '&astra_log\n\n'
    alog += 'AWD       = "%s/"\n'   %awd
    alog += 'AEXT      = "%s/"\n'   %aext
    alog += 'exp_file  = "%s"\n'    %args.exp
    alog += 'equ_file  = "%s"\n'    %args.equ
    alog += 'rev_file  = profiles.dat\n'
    alog += 'tbeg_nml  = %8.4f\n'   %args.tbeg
    alog += 'tend_nml  = %8.4f\n'   %args.tend
    alog += 'TASK      = "%s"\n'    %rtype
    alog += 'machine   = "%s"\n'    %dev_name
    alog += 'debug     = %d\n'      %int(args.debug)
    alog += 'flightsim = %d\n'      %int(args.fs)
    alog += 'resize    = %8.4f\n'   %resize
    if args.tpause is not None:
        alog += 'tpause_nml = %8.4f\n' %args.tpause
    alog += '\n/\n'

    logger.info('Writing %s' %f_log)
    f = open(f_log, 'w')
    f.write(alog)
    f.close()

    logger.debug(alog)

    eqlog = '%s/equ/log/%s' %(awd, args.equ)
    if not os.path.isfile(eqlog):
        input('-----------------------------------\nFile %s missing!\nASTRA will probably crash.\nPress any key to continue at your own risk\n' %eqlog)

    expequ = args.exp + args.equ
    source = '%s/tmp/astra.nml' %awd 
    target = '%s/tmp/%s.nml'    %(awd, expequ) 
    shutil.copy2(source, target)

    cmd = '%s/exe/Build' %awd
    logger.info(cmd)
    os.system(cmd)

    try:
        nc_concat(expequ)
    except:
#        logger.debug(traceback.format_exc())
        pass
