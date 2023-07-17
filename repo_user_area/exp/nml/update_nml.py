src = 'aug34954_lx'

# goblin, tok just as lx

aext_d = { \
    'lx': '/afs/ipp/home/a/astra/ASTRA_LIBRARIES_EXT/rabbit',
    'ldaug': '/shares/departments/AUG/users/git/ASTRA_LIBRARIES_EXT/rabbit',
    'cz': '/compass/home/tardini/astra_ext/rabbit/',
    'ga': '/fusion/projects/codes/astra/ASTRA_LIBRARIES_EXT/rabbit',
    'gway': '/afs/eufus.eu/user/g/g2gtardi/ASTRA_LIBRARIES_EXT/rabbit',
    'iter': '/home/ITER/tarding/ASTRA_LIBRARIES_EXT/rabbit'
}

no_torqjxb_model = ('cz', 'iter', 'gway')

with open(src, 'rb') as f:
    nml = f.read()

for platform in ('tok', 'goblin', 'ldaug', 'ga', 'cz', 'iter', 'gway'):
    dest = 'aug34954_%s' %platform
    nml_dest = nml
    if platform in aext_d.keys():
        nml_dest = nml_dest.replace(aext_d['lx'], aext_d[platform])
    if platform in no_torqjxb_model:
        nml_dest = nml_dest.replace('   torqjxb_model = 3\n', '')

    with open(dest, 'wb') as f:
        f.write(nml_dest)
