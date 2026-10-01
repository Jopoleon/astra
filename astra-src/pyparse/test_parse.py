import os, argparse
from parse_as import LINE2FOR
from equ_parser import EQU_PARSER


def test_parse(statement):

    parse = EQU_PARSER()
    l2f = LINE2FOR(statement, parse)
    print('')
    print(l2f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Write settings & run ASTRA')
    parser.add_argument('-s', '--statement', help='ASTRA statement', required=False)
    args = parser.parse_args()

    if args.statement is None:
        test_parse('HE = MAX(0.,CAR22)+HNGSE')
        test_parse('HE = CAR22+HNGSE')
        test_parse('CMHD1 = THQ99')
    else:
        test_parse(args.statement)
