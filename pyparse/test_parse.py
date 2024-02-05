import os, argparse
import equ_parser
from parse_as import LINE2FOR
from equ_parser import EQU_PARSER

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def test_parse(left, right, f_equ='%s/equ/test' %awd):

    parse = EQU_PARSER(f_equ)
    pack = parse.fml_list, parse.fnc_list, parse.profiles, parse.arr_nam2
    l2f = LINE2FOR(left, right, pack)
    print('')
    print(l2f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Write settings & run ASTRA')
    parser.add_argument('-l', '--left', help='Left-hand side of statement', required=False)
    parser.add_argument('-r', '--right', help='Right-hand side of statement', required=False)
    args = parser.parse_args()

    if args.right is None:
        if args.right is None:
            test_parse('HE', 'MAX(0.,CAR22)+HNGSE')
            test_parse('HE', 'CAR22+HNGSE')
            test_parse('CMHD1', 'THQ99')
        else:
            test_parse('', args.right)
    else:
        if args.right is None:
            test_parse(args.left, '')
        else:
            test_parse(args.left, args.right)
