
id = "VRNT_Gouania_Salaria"

fasta = "alignments/{}_alignment.fa".format(id)
rev_compliments_ids = "reports/reverse_compliments.txt"
blastfiltered = "results/{}_BLAST_filtered.txt".format(id)
outfasta_file = "alignments/{}_alignment_corrected.fa".format(id)

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

with open(fasta, "r") as infasta, open(blastfiltered, "r") as blastrep, open(rev_compliments_ids, "r") as reverseID, open(outfasta_file, "w") as out:

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_IDs = []
    outfastadict = dict()
    blastcoordinates = dict()
    outfasta = str()

    for line in reverseID:
        line = line.rstrip()
        if line.startswith("#"):
            pass
        else:
            line = ">" + line
            reverse_IDs.append(line)

    for id in reverse_IDs:
        for seqID, sequence in read_fasta(infasta):
            shortname = seqID.split(":")[0]
            if shortname == id:
                reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
                id_string = seqID
                seq_string = reverse_complement
                outfastadict[id_string] = seq_string
            else:
                id_string = seqID
                seq_string = sequence
                outfastadict[id_string] = seq_string

    for line in blastrep:
        line = line.rstrip()
        if line.startswith('qseqid'):
            pass
        elif line.split("\t")[8] < line.split("\t")[9]:
            blastcoordinates[str(line.split("\t")[8] + "-" + line.split("\t")[9])] = line.split("\t")[0]
        elif line.split("\t")[8] > line.split("\t")[9]:
            blastcoordinates[str(line.split("\t")[9] + "-" + line.split("\t")[8])] = line.split("\t")[0] # reverse compliment (coordinates changed)

    for seqIDs, sequ in outfastadict.items():
        for blastcoord, queries in blastcoordinates.items():
            if seqIDs.split(":")[1] in blastcoord:
                #print(keys.split(":")[1], queries)
                outfasta += seqIDs + " | query: {}".format(queries) + "\n" + sequ + "\n"

    out.write(outfasta)