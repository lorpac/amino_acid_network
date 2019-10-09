aa_array = [
        ['ALA', 'A'],
        ['ARG', 'R'],
        ['ASN', 'N'],
        ['ASP', 'D'],
        ['CYS', 'C'],
        ['GLU', 'E'],
        ['GLN', 'Q'],
        ['GLY', 'G'],
        ['HIS', 'H'],
        ['ILE', 'I'],
        ['LEU', 'L'],
        ['LYS', 'K'],
        ['MET', 'M'],
        ['PHE', 'F'],
        ['PRO', 'P'],
        ['SER', 'S'],
        ['THR', 'T'],
        ['TRP', 'W'],
        ['TYR', 'Y'],
        ['VAL', 'V']
        ]

three2one_dict = {}
one2three_dict = {}
for aa in aa_array:
    three = aa[0]
    one = aa[1]
    three2one_dict[three] = one
    one2three_dict[one] = three

def three2one(aa):
    try:
        aa_one = three2one_dict[aa]
    except KeyError:
        if len(aa.strip()) == 1:
            aa_one = aa.strip()
    return aa_one

def one2three(aa):
    aa_three = one2three_dict[aa]
    return aa_three