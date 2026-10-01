import os

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
fml_dir = '%s/fml' %awd
fnc_dir = '%s/fnc' %awd

checkeqn = False

NCVA  = 512
NFML  = 500
NSBMX = 20
NRW   = 128

eqn_list = ['NE', 'TE', 'TI', 'CU', 'Equil', 'UPAR', 'F0', 'F1', \
            'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9']

labels_d = { \
    'TE'  : 'Electron temperature', \
    'TI'  : 'Ion temperature', \
    'NE'  : 'Density', \
    'NI'  : 'Ion density', \
    'CU'  : 'Current', \
    'UPAR': 'Parallel velocity'}
for jf in range(10):
    var = 'F%d' %jf
    labels_d[var] = var

short_d = { 'NE': 'N', 'TE': 'E', 'TI': 'I', 'UPAR': 'U', \
    'F0': '0', 'F1': '1', 'F2': '2', 'F3': '3', 'F4': '4', 'F5': '5', \
    'F6': '6', 'F7': '7', 'F8': '8', 'F9': '9' }

coeff_d = { \
    'TE':   ['DE', 'HE', 'XE', 'CE', 'DVE', 'DSE'], \
    'TI':   ['DI', 'HI', 'XI', 'CI', 'DVI', 'DSI'], \
    'NE':   ['DN', 'HN', 'XN', 'CN', 'DVN', 'DSN'], \
    'CU':   ['DC', 'HC', 'XC', 'CC', 'CD'], \
    'UPAR': ['CNPAR', 'XUPAR', 'RUPAR', 'RUPYR', 'RUPFR']}

flux_d = { \
    'TE':   ['PE', 'PET'], \
    'TI':   ['PI', 'PIT'], \
    'NE':   ['SN', 'SNN'], \
    'CU':   ['MV', 'MU'], \
    'UPAR': ['TTRQ']}

bnd_d  = { \
    'TE':   ['TEB', 'QEB', 'QETB'], \
    'TI':   ['TIB', 'QIB', 'QITB'], \
    'NE':   ['NEB', 'QNB', 'QNNB'], \
    'MU':   ['IPL', 'MUB'], \
    'UPAR': ['UPARB', 'TTRQB'] }

for jf in range(10):
    key = 'F%d' %jf
    coeff_d[key] = []
    for lbl in ('DF', 'VF', 'GF', 'DVF', 'DSF'):
        coeff_d[key].append('%s%d' %(lbl, jf))
    flux_d[key] = []
    for lbl in ('SF', 'SFF'):
        flux_d[key].append('%s%d' %(lbl, jf))
    bnd_d[key] = []
    for lbl in ('F', 'QF', 'QFF'):
        bnd_d[key].append('%s%dB' %(lbl, jf))
