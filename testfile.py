def read_fasta(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

with open("alignments/VRNT_Gouania_Salaria_alignment.fa", "r") as input, open('reports/reverse_compliments.txt', "r") as reverseID, open('alignments/VRNT_Gouania_Salaria_alignment_corrected.fa', "w") as out:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    #seq = "TCGGGCCC"
    reverse_IDs = []
    outfasta = str()
    for line in reverseID:
        line = line.rstrip()
        if line.startswith("#"):
            pass
        else:
            line = ">" + line
            reverse_IDs.append(line)
    for id in reverse_IDs:
        for seqID, sequence in read_fasta(input):
            shortname = seqID.split(":")[0]
            if shortname == id:
                reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
                outfasta += seqID + "\n" + reverse_complement + "\n"
            else:
                outfasta += seqID + "\n" + sequence + "\n"
    out.write(outfasta)



