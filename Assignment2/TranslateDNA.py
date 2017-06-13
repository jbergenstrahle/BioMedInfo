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

def reverseComp(seq):
    """Create reverse complement strand from input sequence"""
    compDict = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    comp = ""
    for base in seq:
        comp = compDict[base]+comp
    return comp

def splitSeq(seq):
    """Splits a sequence into codon list"""
    codonsList = {}
    for k in range(0, 3):
        codons = []
        for i in range(k, len(seq)-2, 3):
            codons.append(seq[i:i+3])
        codonsList[k] = codons
    return codonsList

def searchORF(seq):
    """Finds the longest Open Reading Frame (ORF), given a nucleotide sequence.
    ORF is defined as a sequence of codons that can start with all codons (do not need specific start codon)
    and ends with a stop codon (needed). The function returns the protein sequence of the longest ORF"""
    seqForward = seq.upper()

    seqReverse = reverseComp(seq.upper())
    print("Reverse Complementary strand: ", seqReverse)
    #Split sequence into 3-bp lists
    seqForward = splitSeq(seqForward)
    seqReverse = splitSeq(seqReverse)

    ForwardORFStartPos = 0
    ForwardORFEndPos = 0
    ForwardReadingFrame = 0
    ReverseORFStartPos=0
    ReverseORFEndPos = 0
    ReverseReadingFrame = 0
    ORF = []

    #Get the start and stop index
    startAndStopForward, readingFrameORF = getORFIndex(seqForward)
    if startAndStopForward:
        # Look for the longest ORF in forward + reverse sequences
        ForwardORFStartPos, ForwardORFEndPos, ForwardReadingFrame = getStartPosAndLength(startAndStopForward)
        print("Longest ORF Start Position Forward Strand: ",ForwardORFStartPos, " in reading Frame: ",ForwardReadingFrame+1)
    else:
        print("No ORF found in Forward Strand")

    startAndStopReverse, readingFrameORF = getORFIndex(seqReverse)
    if startAndStopReverse:
        ReverseORFStartPos, ReverseORFEndPos, ReverseReadingFrame = getStartPosAndLength(startAndStopReverse)
        print("Longest ORF Start position Reverse Strand: ", ReverseORFStartPos, " in reading Frame: ",ReverseReadingFrame+1)
    else:
        print("No ORF found in Reverse Strand")

    #Evaluate longest between forward and reverse
    if (ForwardORFEndPos - ForwardORFStartPos) > (ReverseORFEndPos - ReverseORFStartPos):
        #Extract ORF from forward strand
        ORF = seqForward[ForwardReadingFrame][ForwardORFStartPos:ForwardORFEndPos]
    else:
        #Extract ORF from reverse strand
        ORF = seqReverse[ReverseReadingFrame][ReverseORFStartPos:ReverseORFEndPos]

    if len(ORF)>0:
        print("Longest ORF found: ", ORF)
        proteinSeq = translateProtein(ORF)
        return proteinSeq
    else:
        print("No protein sequence present")

def getStartPosAndLength(seqIdx):
    """Finds the longest ORF from a list of start and stop indicies
     ([start1_RF1, start2_RF1] , [stop1_RF1, stop2_FR1]), [start1_RF1, start2_RF2 .... etc)
     returns a list with start and stop and readingframe for the longest ORF"""
    #print("seqIdx", seqIdx)
    lenORF = 0
    readingFrameOfLongestORF = 0
    longestORFStartPos = 0
    longestORFEndPos = 0

   # print("length: ", len(seqIdx))
    for idx, RF in enumerate(seqIdx): #adds index to each value -> enumerate (idx gets enumerate value, RF gets each seqIdx value)

        for ORF in zip(*RF): #* will "unlist" seqIdx list for the zip function zip(*seqIdx)=zip(L1,L2)
            #print("ORF:", ORF)
            currLenORF = ORF[1]-ORF[0]
            if currLenORF > lenORF:
            #Which correponds to this ORF:
                lenORF = currLenORF
                longestORFStartPos = ORF[0]
                longestORFEndPos = ORF[1]
                readingFrameOfLongestORF = idx

    return longestORFStartPos, longestORFEndPos, readingFrameOfLongestORF

def getORFIndex(codonList):
    """Finds all ORF given a list of codons"""
    startAndStopIdx = [[],[],[]]
#Loop through each forward and reverse list
    for n in codonList:
        codons = codonList[n]
        #print("codons: ", codons)
        readingFrameOfORFs = []

        stopIdx = []
        startIdx = []
        idx = 0
        #start from beginning (if stop is first codon, I get first stop at same place, i.e. 0 length ORF)
        startIdx.append(idx)

        for codon in codons:
            if codon in stopCodons:
                stopIdx.append(idx)
                #print("stopIdx",stopIdx)
                startIdx.append(idx+1) #next start codon is next after the last stop codon
            idx += 1

        startIdx.pop() #remove last entry (out of sequence)
        if len(stopIdx) == 0:
           print("no stop codon found")
           #stopIdx = 0
        elif not startIdx:
           print("no stop codon found after first codon")
           #startIdx = 0
        else:
            print("ORF found in reading frame ", n+1)
            readingFrameOfORFs.append(n)
            startAndStopIdx[n].append(startIdx)
            startAndStopIdx[n].append(stopIdx)

    return startAndStopIdx, readingFrameOfORFs

def translateProtein(ORF):
    """Translate nucleotide sequence to protein sequence"""
    protseq = ""
    for codon in ORF:
        protseq = aminoAcidMap[codon] + protseq
    return protseq

#variables
stopCodons = ["TAA", "TGA", "TAG"]

aminoAcidMap = {'TTT': 'F',
                'TTC': 'F',
                'TTA': 'L',
                'TTG': 'L',

                'TCT': 'S',
                'TCC': 'S',
                'TCA': 'S',
                'TCG': 'S',

                'TAT': 'Y',
                'TAC': 'Y',
                'TAA': '*',
                'TAG': '*',

                'TGT': 'C',
                'TGC': 'C',
                'TGA': '*',
                'TGG': 'W',

                'CTT': 'L',
                'CTC': 'L',
                'CTA': 'L',
                'CTG': 'L',

                'CCT': 'P',
                'CCC': 'P',
                'CCA': 'P',
                'CCG': 'P',

                'CAT': 'H',
                'CAC': 'H',
                'CAA': 'Q',
                'CAG': 'Q',

                'CGT': 'R',
                'CGC': 'R',
                'CGA': 'R',
                'CGG': 'R',

                'ATT': 'I',
                'ATC': 'I',
                'ATA': 'I',
                'ATG': 'M',

                'ACT': 'T',
                'ACC': 'T',
                'ACA': 'T',
                'ACG': 'T',

                'AAT': 'N',
                'AAC': 'N',
                'AAA': 'K',
                'AAG': 'K',

                'AGT': 'S',
                'AGC': 'S',
                'AGA': 'R',
                'AGG': 'R',

                'GTT': 'V',
                'GTC': 'V',
                'GTA': 'V',
                'GTG': 'V',

                'GCT': 'A',
                'GCC': 'A',
                'GCA': 'A',
                'GCG': 'A',

                'GAT': 'D',
                'GAC': 'D',
                'GAA': 'E',
                'GAG': 'E',

                'GGT': 'G',
                'GGC': 'G',
                'GGA': 'G',
                'GGG': 'G'}

testSeq = "ATTCGCGCGCGTAAGCGCCCCCCCCCTAAGGCGGTTAAAATTTGGCCCGGAGCGAGCGTGGTCGGTTT"
testSeq2 = "GCTCGTAATGATTGTCAAGAAGGTCATATTTTAAAAATGTTTCCTTCTACTTGGTATGTT"
testSeq3 = "gctcgtaatgattgtcaagaaggtcatattttaaaaatgtttccttctacttggtatgtt"
testSeq4 = "TAACAGTAG"
testSeq5 = "A"
testSeq6 = "NNNNNNNNNNNNNNNNN"
testSeq7 = "TAACGTAGNNNCCGN"

parser = argparse.ArgumentParser() #automatic adds -h option
parser.add_argument("inFasta") #positionellt argument, den kommer först
parser.add_argument("-v", "--version", action="version", version="1.0")
args = parser.parse_args()


fasta = args.inFasta
ids, seqs = readFastq(fasta)
protSeq = searchORF(seqs[0])
#protSeq = searchORF(testSeq7)
print("Protein Sequence: ", protSeq)