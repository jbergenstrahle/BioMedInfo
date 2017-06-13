import argparse

def readFastq(fasta):
    """Parse Fasta files"""
    ids = []
    seqs = []
    with open(fasta) as fa:
        while True: #loop until end of file
            line = fa.readline().rstrip() #rstrip to strip whitespace
            if len(line) == 0:
                break
            if line[0] == '>':
                ids.append(line[1:])
                seqs.append('')
            else:
                seqs[-1] += line.upper()
    return ids, seqs

parser = argparse.ArgumentParser()
parser.add_argument("-in", required=True)
parser.add_argument("-motif", required=True)
args = parser.parse_args()

