import os
import subprocess


def create_par(gtf, gen, tmp, n, l, wd):
    par_path = os.path.join(wd, "simulation.par")
    with open(par_path, "w") as par:
        par.write(f"REF_FILE_NAME\t{gtf}\n")
        par.write(f"GEN_DIR\t{gen}\n")
        par.write(f"TMP_DIR\t{tmp}\n")
        par.write("POLYA_SCALE\tNaN\n")  # TODO: add to cli
        par.write("POLYA_SHAPE\tNaN\n")  # TODO: add to cli
        par.write(f"READ_NUMBER\t{n}\n")
        par.write(f"READ_LENGTH\t{l}\n")
        par.write("PAIRED_END\tYES\n")  # TODO: add to cli
        par.write("FASTA\tYES\n")
        err_model = "76" if l >= 76 else "35"
        par.write(f"ERR_FILE\t{err_model}\n")


def run_flux(wd, threads):
    par_path = os.path.join(wd, "simulation.par")
    fluxlog_path = os.path.join(wd, "flux.log")
    fluxout_path = os.path.join(wd, "flux.out")
    flux_cmd = [
        "flux-simulator",
        "-t",
        "simulator",
        "-x",
        "-l",
        "-s",
        "-p",
        par_path,
        "--threads",
        f"{threads}",
    ]
    flux = subprocess.run(
        flux_cmd,
        stdout=open(fluxout_path, "w"),
        stderr=open(fluxlog_path, "w"),
    )
    return flux.returncode, fluxlog_path


def run_flux_x2(wd1, wd2, threads):
    threads = max(1, threads / 2)

    par1_path = os.path.join(wd1, "simulation.par")
    fluxlog1_path = os.path.join(wd1, "flux.log")
    fluxout1_path = os.path.join(wd1, "flux.out")
    par2_path = os.path.join(wd2, "simulation.par")
    fluxlog2_path = os.path.join(wd2, "flux.log")
    fluxout2_path = os.path.join(wd2, "flux.out")

    flux1_cmd = [
        "flux-simulator",
        "-t",
        "simulator",
        "-x",
        "-l",
        "-s",
        "-p",
        par1_path,
        "--threads",
        f"{threads}",
    ]
    flux2_cmd = [
        "flux-simulator",
        "-t",
        "simulator",
        "-x",
        "-l",
        "-s",
        "-p",
        par2_path,
        "--threads",
        f"{threads}",
    ]

    flux1 = subprocess.Popen(
        flux1_cmd,
        stdout=open(fluxout1_path, "w"),
        stderr=open(fluxlog1_path, "w"),
    )
    flux2 = subprocess.Popen(
        flux2_cmd,
        stdout=open(fluxout2_path, "w"),
        stderr=open(fluxlog2_path, "w"),
    )
    flux1.wait()
    flux2.wait()

    return flux1.returncode, fluxlog1_path, flux2.returncode, fluxlog2_path
