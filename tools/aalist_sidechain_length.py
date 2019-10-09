import amino_acids_conversion as aac
import operator

dict_size = {
        'A': 1.5,
        'C': 2.8,
        'D': 3.7,
        'E': 4.6,
        'F': 5.0,
        'G': 0.0,
        'H': 4.6,
        'I': 3.9,
        'K': 5.9,
        'L': 3.7,
        'M': 4.7,
        'N': 3.5,
        'P': 2.4,
        'Q': 4.6,
        'R': 6.5,
        'S': 2.4,
        'T': 2.5,
        'V': 2.5,
        'W': 6.6,
        'Y': 6.4}

list_size = sorted(dict_size.items(), key=operator.itemgetter(1))

with open("aminoacids_sidechain.txt", "w") as f:
        for aa, _ in list_size:
                aa = aac.one2three(aa)
                f.write(aa)
                f.write("\n")

