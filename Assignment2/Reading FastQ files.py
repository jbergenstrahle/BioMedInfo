def readFastq(fasta):
    """Parse FastQ files"""
    ids = []
    seqs = []
    dess = []
    quals = []
    with open(fasta) as fa:
        while True: #loop until end of file
            id = fa.readline().rstrip() #rstrip to strip whitespace
            seq = fa.readline().rstrip()
            des = fa.readline().rstrip()
            qual = fa.readline().rstrip()
            if len(id) == 0:
                break
            ids.append(id)
            seqs.append(seq)
            dess.append(des)
            quals.append(qual)
    return ids, seqs, dess, quals

ids, seqs, dess, quals = readFastq("fasta.fasta")

print(len(ids))
for line in ids:
    print(line)



