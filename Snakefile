import pandas as pd
import numpy as np
import os
jn = os.path.join

## Config
configfile: "config.yaml"

### Usage

rule all:
    input:
         "alignments/{}_alignment_corrected_aligned.fa".format(config["blast"]["query_name"])

rule create_local_blastdb:
    input:
        config["blast"]["fasta_db"]
    output:
        "blast_db/checkpoint_blast_db.txt"
    params:
        db_name = "blast_db/{}".format(config["blast"]["db_name"])
    shell:
        """
        makeblastdb -in {input} -out {params.db_name} -dbtype nucl
        touch {output}
        """

rule blastn:
    input:
        "blast_db/checkpoint_blast_db.txt",
        query = config["blast"]["query_fasta"]
    output:
        "reports/{}_BLAST.txt".format(config["blast"]["query_name"])
    params:
        tmp_out = "blast_tmp",
        db_name = "blast_db/{}".format(config["blast"]["db_name"])
    shell:
        """
        blastn -query {input.query} -db {params.db_name} -out {params.tmp_out} -outfmt 6
        echo $'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' | cat - {params.tmp_out} > {output}
        rm {params.tmp_out}
        """

rule filter_blast_report:
    input:
        "reports/{}_BLAST.txt".format(config["blast"]["query_name"])
    output:
        "reports/{}_BLAST_filtered.txt".format(config["blast"]["query_name"])
    params:
        e_value = config["filter_blast"]["evalue"],
        identity = config["filter_blast"]["identity"]
    run:
        fh = pd.read_csv(input[0], sep = "\t", index_col=0)
        fh = fh.loc[(fh['pident'] > int(params.identity)) & (fh['evalue'] < float(params.e_value))]
        fh.to_csv(output[0], sep = "\t")

rule bedtools_prep:
    input:
        "reports/{}_BLAST_filtered.txt".format(config["blast"]["query_name"])
    output:
        bed = "reports/{}_BLAST_filtered.bed".format(config["blast"]["query_name"]),
        rev_compliments_ids = "reports/{}_reverse_compliments.txt".format(config["blast"]["query_name"])
    run:
        with open(input[0], "r") as blastreport, open(output.bed, "w") as bedinput, open(output.rev_compliments_ids, "w") as revcom_IDs:
            bed = str()
            revcom = "## This file contains all sequence IDs that need to be reverse complemented before an alignment step.\n"
            for line in blastreport:
                if line.startswith('qseqid'):
                    pass
                elif int(line.split("\t")[8]) < int(line.split("\t")[9]):
                    bed += line.split("\t")[1] + "\t" + line.split("\t")[8] + "\t" + line.split("\t")[9] + '\n' # normal the start is smaller than the end coordinate
                elif int(line.split("\t")[8]) > int(line.split("\t")[9]):
                    revcom += line.split("\t")[1] + '\n'
                    bed += line.split("\t")[1] + "\t" + line.split("\t")[9] + "\t" + line.split("\t")[8] + '\n' # reverse complement needed the end is smaller than the start coordinate
            bedinput.write(bed)
            revcom_IDs.write(revcom)

rule bedtools_getfasta:
    input:
        bed = "reports/{}_BLAST_filtered.bed".format(config["blast"]["query_name"]),
        fasta_db = config["blast"]["fasta_db"]
    output:
        fasta = "alignments/{}_alignment.fa".format(config["blast"]["query_name"])
    conda:
        "envs/bedtools.yml"
    shell:
        "bedtools getfasta -fi {input.fasta_db} -bed {input.bed} -fo {output.fasta}"


rule reverse_complement:
    input:
        fasta = "alignments/{}_alignment.fa".format(config["blast"]["query_name"]),
        rev_compliments_ids = "reports/{}_reverse_compliments.txt".format(config["blast"]["query_name"]),
        blastfiltered = "reports/{}_BLAST_filtered.txt".format(config["blast"]["query_name"])
    output:
        outfasta = "alignments/{}_alignment_corrected.fa".format(config["blast"]["query_name"])
    run:
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

        with open(input.fasta, "r") as infasta, open(input.blastfiltered, "r") as blastrep, open(input.rev_compliments_ids, "r") as reverseID, open(output.outfasta, "w") as out:
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
                    blastcoordinates.setdefault(str(line.split("\t")[8] + "-" + line.split("\t")[9]), []).append(line.split("\t")[0])
                elif line.split("\t")[8] > line.split("\t")[9]:
                    blastcoordinates.setdefault(str(line.split("\t")[9] + "-" + line.split("\t")[8]), []).append(line.split("\t")[0]) # reverse compliment (coordinates changed)

            for seqIDs, sequ in outfastadict.items():
                for blastcoord, queries in blastcoordinates.items():
                    if seqIDs.split(":")[1] in blastcoord:
                        outfasta += seqIDs + " | blast-query: {}".format(queries) + "\n" + sequ + "\n"

            out.write(outfasta)

rule muscle_align:
    input:
        "alignments/{}_alignment_corrected.fa".format(config["blast"]["query_name"])
    output:
        "alignments/{}_alignment_corrected_aligned.fa".format(config["blast"]["query_name"])
    conda:
        "envs/muscle.yml"
    params:
        iterations = config["muscle_alignment"]["iterations"]
    shell:
         "muscle -in {input} -out {output} -maxiters {params.iterations}"
