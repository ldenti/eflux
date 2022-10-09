import os
from Bio import SeqIO
import gffutils


def split_fa(fapath, odir):
    for record in SeqIO.parse(fapath, "fasta"):
        outfa = open(os.path.join(odir, "{}.fa".format(record.id)), "w")
        SeqIO.write(record, outfa, "fasta")
        outfa.close()


def open_gtf(gtf_path):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            dbfn="{}.db".format(gtf_path),
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return gtf


def split_fq(reads_path, ofq1_path, ofq2_path):
    ofq1 = open(ofq1_path, "w")
    ofq2 = open(ofq2_path, "w")
    records1 = []
    records2 = []
    for record in SeqIO.parse(reads_path, "fastq"):
        if record.id[-1] == "1":
            records1.append(record)
        else:
            records2.append(record)
        if len(records1) > 5000:
            SeqIO.write(records1, ofq1, "fastq")
            records1 = []
        if len(records2) > 5000:
            SeqIO.write(records2, ofq2, "fastq")
            records2 = []
    if len(records1) > 0:
        SeqIO.write(records1, ofq1, "fastq")
        records1 = []
    if len(records2) > 0:
        SeqIO.write(records2, ofq2, "fastq")
        records2 = []
    ofq1.close()
    ofq2.close()
