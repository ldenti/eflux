import argparse


def get_args():
    parser = argparse.ArgumentParser(description="ASFlux")
    parser.add_argument(
        "FA",
        type=str,
        help="Path to FASTA reference",
    )
    parser.add_argument(
        "GTF",
        type=str,
        help="Path to GTF annotation",
    )

    parser.add_argument(
        "-e",
        "--events",
        dest="events",
        default="SE,SS,MX,RI",
        type=str,
        help="Comma separated list of events (default: SE,SS,MX,RI)",
    )
    parser.add_argument(
        "-l",
        "--length",
        dest="l",
        default=100,
        type=int,
        help="Read length (default: 100)",
    )
    parser.add_argument(
        "-n",
        "--nreads",
        dest="n",
        default=1000,
        type=int,
        help="Number of reads (default: 1000)",
    )
    parser.add_argument(
        "-p",
        "--prob",
        dest="P",
        default=100,
        type=int,
        help="Probability to select a gene (default: 100)",
    )

    parser.add_argument(
        "--simulation",
        dest="simulation",
        default=False,
        action="store_true",
        help="Run only simulation",
    )

    # parser.add_argument(
    #     "-s",
    #     "--seed",
    #     dest="seed",
    #     default=23,
    #     type=int,
    #     help="Seed for random generator (default: 23)",  # TODO: use the current time
    # )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        default=1,
        type=int,
        help="Number of threads (default: 1)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Verbose mode",
    )
    parser.add_argument(
        "--version",
        dest="version",
        action="store_true",
        default=False,
        help="Print version",
    )

    # Directories
    parser.add_argument(
        "--out",
        dest="wd",
        default="asflux-wd",
        type=str,
        help="Output directory (default: asflux-wd)",
    )
    parser.add_argument(
        "--tmp",
        dest="tmpd",
        default="asflux-tmp",
        type=str,
        help="Temporary directory (default: asflux-tmp)",
    )

    return parser.parse_args()