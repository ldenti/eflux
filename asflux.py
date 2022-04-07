import sys, os
import argparse
import logging
from datetime import datetime
import subprocess

from Bio import SeqIO


def time():
    return datetime.now().strftime("%B %d %Y - %H:%M:%S")


def split_fa(fapath, odir):
    for record in SeqIO.parse(fapath, "fasta"):
        outfa = open(os.path.join(odir, "{}.fa".format(record.id)), "w")
        SeqIO.write(record, outfa, "fasta")
        outfa.close()


def main():
    parser = argparse.ArgumentParser(description="asflux")

    parser.add_argument("FA", type=str, help="Path to reference (FASTA)")
    parser.add_argument("GTF", type=str, help="Path to annotation (GTF)")
    parser.add_argument(
        "--wd",
        dest="wd",
        type=str,
        default=".",
        help="Path to working directory (default:.)",
    )
    parser.add_argument(
        "-l",
        "--length",
        dest="l",
        type=int,
        default=100,
        help="Read length (default: 100)",
    )
    parser.add_argument(
        "-r",
        "--reads",
        dest="r",
        type=int,
        default=5000000,
        help="Number of reads (default: 5000000)",
    )
    parser.add_argument(
        "--tmp",
        dest="tmpd",
        type=str,
        default=".",
        help="Path to temporary directory (default:.)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=2,
        help="Number of threads (default: 2)",
    )
    parser.add_argument("--debug", action="store_true", help="Activate debug mode")

    args = parser.parse_args()

    FORMAT = f"[{time()}] %(message)s"
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.DEBUG if args.debug else logging.INFO,
        format=FORMAT,
    )

    print("")
    logging.info("Setting up..")
    args.wd = os.path.abspath(args.wd)
    args.GTF = os.path.abspath(args.GTF)
    if args.wd != "." and os.path.isdir(args.wd) and len(os.listdir(args.wd)) != 0:
        logging.error(f"Working directory {args.wd} is not empty. Halting.\n")
        sys.exit(1)
    os.makedirs(args.wd, exist_ok=True)
    par_path = os.path.join(args.wd, "simulation.par")
    refs_path = os.path.join(args.wd, "references")
    os.makedirs(refs_path, exist_ok=True)
    fluxlog_path = os.path.join(args.wd, "flux.log")
    ofq1_path = os.path.join(args.wd, "sample_1.fq")
    ofq2_path = os.path.join(args.wd, "sample_2.fq")

    logging.info("Splitting the reference..")
    split_fa(args.FA, refs_path)

    logging.info("Populating PAR file..")
    with open(par_path, "w") as par:
        par.write(f"REF_FILE_NAME\t{args.GTF}\n")
        par.write(f"GEN_DIR\t{refs_path}\n")
        par.write(f"TMP_DIR\t{args.tmpd}\n")
        par.write("POLYA_SCALE\tNaN\n")
        par.write("POLYA_SHAPE\tNaN\n")
        par.write(f"READ_NUMBER\t{args.r}\n")
        par.write(f"READ_LENGTH\t{args.l}\n")
        par.write("PAIRED_END\tYES\n")
        par.write("FASTA\tYES\n")
        err_model = "76" if args.l >= 76 else "35"
        par.write(f"ERR_FILE\t{err_model}\n")

    logging.info("Running flux simulator..")
    flux = subprocess.run(
        [
            "flux-simulator",
            "-t",
            "simulator",
            "-x",
            "-l",
            "-s",
            "-p",
            par_path,
            "--threads",
            str(args.threads),
        ],
        stdout=open("/dev/null", "w"),
        stderr=open(fluxlog_path, "w"),
    )
    if flux.returncode > 0:
        logging.error(
            f"Flux crashed (return code {flux.returncode}). See {fluxlog_path} for more info. Halting..\n"
        )
        sys.exit(1)

    logging.info("Splitting reads..")
    reads_path = os.path.join(args.wd, "simulation.fastq")
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
    logging.info(f"Samples: {ofq1_path} and {ofq2_path}.")
    logging.info("Done.\n")


if __name__ == "__main__":
    main()
