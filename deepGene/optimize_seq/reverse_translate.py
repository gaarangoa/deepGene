import json

d_aa_nt = { #inverse translation (what an amino acid can be coded into a nt sequence codon)
    'M':[{0:{'A'}, 1:{'T'}, 2:{'G'}}]
    , 'F':[{0:{'T'}, 1:{'T'}, 2:{'T','C'}}]
    , 'L':[{0:{'T'}, 1:{'T'}, 2:{'A', 'G'}}, {0:{'C'}, 1:{'T'}, 2:{'A', 'C', 'G','T'}}]
    , 'I':[{0:{'A'}, 1:{'T'}, 2:{'A', 'T', 'C'}}]
    , 'V':[{0:{'G'}, 1:{'T'}, 2:{'A', 'C', 'G','T'}}]
    , 'S':[{0:{'T'}, 1:{'C'}, 2:{'A', 'C', 'G','T'}}, {0:{'A'}, 1:{'G'}, 2:{'T','C'}}]
    , 'P':[{0:{'C'}, 1:{'C'}, 2:{'A', 'C', 'G','T'}}]
    , 'T':[{0:{'A'}, 1:{'C'}, 2:{'A', 'C', 'G','T'}}]
    , 'A':[{0:{'G'}, 1:{'C'}, 2:{'A', 'C', 'G','T'}}]
    , 'Y':[{0:{'T'}, 1:{'A'}, 2:{'T', 'C'}}]
    , 'H':[{0:{'C'}, 1:{'A'}, 2:{'T', 'C'}}]
    , 'Q':[{0:{'C'}, 1:{'A'}, 2:{'A', 'G'}}]
    , 'N':[{0:{'A'}, 1:{'A'}, 2:{'T', 'C'}}]
    , 'K':[{0:{'A'}, 1:{'A'}, 2:{'A', 'G'}}]
    , 'D':[{0:{'G'}, 1:{'A'}, 2:{'T', 'C'}}]
    , 'E':[{0:{'G'}, 1:{'A'}, 2:{'A', 'G'}}]
    , 'C':[{0:{'T'}, 1:{'G'}, 2:{'T', 'C'}}]
    , 'R':[{0:{'C'}, 1:{'G'}, 2:{'A', 'C', 'G','T'}}, {0:{'A'}, 1:{'G'}, 2:{'A', 'G'}}]
    , 'G':[{0:{'G'}, 1:{'G'}, 2:{'A', 'C', 'G','T'}}]
    , 'W':[{0:{'T'}, 1:{'G'}, 2:{'G'}}]
}

@staticmethod
def dump_dict(file_name):
    with open(file_name, 'w') as f:
        json.dump(d_aa_nt, f)

def get_dict():
    return d_aa_nt
