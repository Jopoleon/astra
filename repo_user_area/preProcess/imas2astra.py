import os, argparse
from trview import plasma_state

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def imas2astra(shot, idsRun):
    plasma = plasma_state.PLASMA_STATE(shot)
    plasma.fromIMAS(idsRun)
    plasma.toASTRA(awd=awd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reading IMAS output and writing ASTRA input')    
    parser.add_argument('-s', '--shot', type=int, help='Shot number', required=False, default=38384)
    parser.add_argument('-r', '--ids_run', type=int, help='#IDS run', required=False, default=6)
    args = parser.parse_args()
    imas2astra(args.shot, args.ids_run)
