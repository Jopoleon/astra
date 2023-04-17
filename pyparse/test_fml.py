import equ_parser
from parse_as import *

awd='/afs/ipp/home/g/git/a82'
f_equ = '%s/tmp/model.tmp' %awd
parse = equ_parser.EQU_PARSER(f_equ)

pack = parse.fml_list, parse.fnc_list, parse.profiles, parse.arr_nam2

#l2f = LINE2FOR('CC', 'cnhr + nuee', pack)

l2f = LINE2FOR('HE', 'MAX(0.,CAR22)+HNGSE', pack)
print('')
print(l2f.fcode)

l2f = LINE2FOR('HE', 'CAR22+HNGSE', pack)
print('')
print(l2f.fcode)

l2f = LINE2FOR('CMHD1', 'THQ99', pack)
print('')
print(l2f.fcode)
