import os

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def IncludedFML(fml_in):

    fml_out = []
    for formula in fml_in:
        ffml = '%s/fml/%s' %(awd, formula.lower())
        with open(ffml, 'r') as f:
            lines = f.readlines()

        for line in lines:
            line = line.lower().strip()
            if 'include' in line and line[0] != '!':
                fml = line.split('fml/')[1].split('\'')[0]
                fml_out.append(fml)
    return fml_out


class FML:


    def __init__(self, pieces, fml_list):

        self.allFML(pieces, fml_list)
        self.textBlock()


    def allFML(self, pieces, fml_list):

        self.fml_all = []
        fmls = [var.lower() for var in pieces if var.upper() in fml_list]
        self.rec_fml(fmls)
        self.fml_all = self.fml_all[::-1] + fmls
        self.fml_unique = []
        for x in self.fml_all:
            if x not in self.fml_unique:
                self.fml_unique.append(x)


    def rec_fml(self, fml_in):

        fmll = IncludedFML(fml_in)
        if fmll:
            self.fml_all += fmll
            self.rec_fml(fmll)


    def textBlock(self):

        self.txt = ''
        for formula in self.fml_unique:
            ffml = '%s/fml/%s' %(awd, formula)
            with open(ffml, 'r') as f:
                lines = f.readlines()
            for line in lines:
                line = line.upper().strip()
                if line:
                    if line[0] != '!' and 'INCLUDE' not in line:
                        self.txt += line + '\n'


if __name__ == '__main__':

    import argparse
    from parse_as import rec_split
    from equ_parser import EQU_PARSER

    parser = argparse.ArgumentParser(description='Write settings & run ASTRA')
    parser.add_argument('-s', '--statement', help='ASTRA statement', required=False)
    parser.add_argument('-f', '--equ', help='ASTRA equ file', required=False, default='%s/equ/test' %awd)
    args = parser.parse_args()

    parse = EQU_PARSER(args.equ)
    if args.statement is None:
        pieces = rec_split('HE = MAX(0.,CAR22)+HNGSE+ CNSA + 3.*PEDT')
    else:
        pieces = rec_split(args.statement)

    fml = FML(pieces, parse.fml_list)

    print(fml.txt)
