import sys
import os
import logging

from cli import get_args
from utils import split_fa, open_gtf, split_fq
from events import run_suppa, get_events, select_events, split_annotation
from reads import create_par, run_flux

FORMAT = "[%(asctime)s] %(message)s"
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)

VERSION = "0.0.2"


def main():
    if "--version" in sys.argv:
        print(f"eflux v{VERSION}")
        sys.exit(0)

    args = get_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    args.wd = os.path.join(os.getcwd(), args.wd)
    args.tmpd = os.path.join(os.getcwd(), args.tmpd)

    try:
        os.makedirs(args.wd)
    except FileExistsError:
        logging.critical("Output folder already exits.")
        logging.critical("Halting..\n")
        sys.exit(1)
    if not os.path.isdir(args.tmpd):
        os.makedirs(args.tmpd)

    gtf1_path = args.GTF
    gtf2_path = ""
    if not args.simulation:
        logging.info("Running SUPPA2 to generate events..")
        suppa2_wd = os.path.join(args.wd, "suppa2")
        os.makedirs(suppa2_wd)
        retcode, suppa2log_path = run_suppa(args.GTF, suppa2_wd, args.events.split(","))
        if retcode != 0:
            logging.critical(f"SUPPA2 did not run succesfully (return code {retcode}).")
            logging.critical(f"See {suppa2log_path} for more details.")
            logging.critical("Halting..\n")
            sys.exit(1)
        logging.critical("SUPPA2 ran succesfully.")
        logging.info("Getting events from SUPPA2 output..")
        events = get_events(suppa2_wd)
        logging.info("Selecting events..")
        events = select_events(events, args.P)

        logging.info("Opening input annotation..")
        gtf = open_gtf(args.GTF)

        logging.info("Splitting annotations..")
        gtf1_path, gtf2_path = split_annotation(gtf, args.wd, events)

    chroms_path = args.chroms
    if chroms_path == None:
        chroms_path = os.path.join(args.wd, "references")
        logging.info(f"Splitting reference to {chroms_path}..")
        os.makedirs(chroms_path, exist_ok=True)
        split_fa(args.FA, chroms_path)
    else:
        logging.info(f"Using {chroms_path}..")

    flux_wd = os.path.join(args.wd, "flux")
    os.makedirs(flux_wd, exist_ok=True)

    logging.info("Populating PAR file..")
    create_par(gtf1_path, chroms_path, args.tmpd, args.n, args.l, flux_wd)

    logging.info("Running flux simulator..")
    retcode, fluxlog_path = run_flux(flux_wd, args.threads)
    if retcode > 0:
        logging.critical(f"Flux did not run succesfully (return code {retcode}).")
        logging.critical(f"See {fluxlog_path} for more details.")
        logging.critical("Halting..\n")
        sys.exit(1)

    reads_path = os.path.join(flux_wd, "simulation.fastq")
    ofq1_path = os.path.join(args.wd, "sample_1.fq")
    ofq2_path = os.path.join(args.wd, "sample_2.fq")

    logging.info("Splitting reads..")
    split_fq(reads_path, ofq1_path, ofq2_path)
    logging.info(f"Samples: {ofq1_path} and {ofq2_path}.")

    logging.info("Done.\n")


if __name__ == "__main__":
    main()
