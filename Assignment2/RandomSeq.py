import random

seq=""

while True:
    try:
        seqLength = int(input("Enter Sequence Length: "))
        break
    except ValueError:
        print ("That is not a number... Try again")


for i in range(seqLength):
    seq += random.choice("ACGT")

#create fasta
fasta = ">randomID\n{}".format(seq) #ersätt { } med seq

print(fasta)


