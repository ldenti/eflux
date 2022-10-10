import sys
import os
import logging
from multiprocessing import Pool  # , set_start_method

from cli import get_args
from utils import split_fa, open_gtf, split_fq
from events import run_suppa, get_events, select_events, split_annotation
from reads import create_par, run_flux, run_flux_x2

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

    c1gtf_path = args.GTF
    c2gtf_path = ""
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
        c1gtf_path, c2gtf_path = split_annotation(gtf, args.wd, events)

    chroms_path = args.chroms
    if chroms_path is None:
        chroms_path = os.path.join(args.wd, "references")
        logging.info(f"Splitting reference to {chroms_path}..")
        os.makedirs(chroms_path, exist_ok=True)
        split_fa(args.FA, chroms_path)
    else:
        logging.info(f"Using chromosomes in {chroms_path}..")

    for replicate in range(1, args.replicates + 1):
        flux1_wd = os.path.join(args.wd, f"flux-{replicate}-C1")
        flux1_tmp = os.path.join(args.tmpd, f"flux-{replicate}-C1")
        os.makedirs(flux1_wd)
        os.makedirs(flux1_tmp)
        logging.info(f"Replicate {replicate}: Populating PAR file..")
        create_par(c1gtf_path, chroms_path, flux1_tmp, args.n, args.l, flux1_wd)
        if c2gtf_path != "":
            flux2_wd = os.path.join(args.wd, f"flux-{replicate}-C2")
            flux2_tmp = os.path.join(args.tmpd, f"flux-{replicate}-C2")
            os.makedirs(flux2_wd)
            os.makedirs(flux2_tmp)
            create_par(c2gtf_path, chroms_path, flux2_tmp, args.n, args.l, flux2_wd)

        logging.info(f"Replicate {replicate}: Running flux simulator..")
        retcode1, fluxlog1_path, retcode2, fluxlog2_path = 0, "", 0, ""
        if c2gtf_path == "":
            retcode1, fluxlog1_path = run_flux(flux1_wd, args.threads)
        else:
            retcode1, fluxlog1_path, retcode2, fluxlog2_path = run_flux_x2(
                flux1_wd, flux2_wd, args.threads
            )

        if retcode1 > 0:
            logging.critical(f"Flux did not run succesfully (return code {retcode1}).")
            logging.critical(f"See {fluxlog1_path} for more details.")
            logging.critical("Halting..\n")
            sys.exit(1)
        if retcode2 > 0:
            logging.critical(f"Flux did not run succesfully (return code {retcode2}).")
            logging.critical(f"See {fluxlog2_path} for more details.")
            logging.critical("Halting..\n")
            sys.exit(1)

        logging.info("Splitting reads..")
        with Pool(processes=2) as pool:
            reads1_path = os.path.join(flux1_wd, "simulation.fastq")
            ofq11_path = os.path.join(args.wd, f"R{replicate}-C1_1.fq")
            ofq12_path = os.path.join(args.wd, f"R{replicate}-C1_2.fq")
            reads2_path = os.path.join(flux2_wd, "simulation.fastq")
            ofq21_path = os.path.join(args.wd, f"R{replicate}-C2_1.fq")
            ofq22_path = os.path.join(args.wd, f"R{replicate}-C2_2.fq")
            _ = pool.starmap(
                split_fq,
                [
                    (reads1_path, ofq11_path, ofq12_path),
                    (reads2_path, ofq21_path, ofq22_path),
                ],
            )
            logging.info(
                f"Replicate/Condition {replicate}-1: {ofq11_path} and {ofq12_path}."
            )
            logging.info(
                f"Replicate-Condition {replicate}-2: {ofq21_path} and {ofq22_path}."
            )

    logging.info("Done.\n")


if __name__ == "__main__":
    main()
