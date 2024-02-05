import equ_parser
from parse_as import *

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def test_parse(left, right, f_equ='%s/equ/test' %awd):

    parse = equ_parser.EQU_PARSER(f_equ)
    pack = parse.fml_list, parse.fnc_list, parse.profiles, parse.arr_nam2
    l2f = LINE2FOR(left, right, pack)
    print('')
    print(l2f)


if __name__ == '__main__':

    test_parse('HE', 'MAX(0.,CAR22)+HNGSE')
    test_parse('HE', 'CAR22+HNGSE')
    test_parse('CMHD1', 'THQ99')
