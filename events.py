import os
import glob
import random
import subprocess


def run_suppa(GTF, WD, events):
    suppa2log_path = os.path.join(WD, "suppa2.log")
    suppa2_cmd = [
        "suppa.py",
        "generateEvents",
        "-i",
        GTF,
        "-f",
        "ioe",
        "-o",
        WD + "/",
        "-e",
    ] + events
    p = subprocess.run(suppa2_cmd, stderr=open(suppa2log_path, "w"))
    return p.returncode, suppa2log_path


def get_events(WD):
    events = {}
    for ioe in glob.glob(os.path.join(WD, f"_*_strict.ioe")):
        for line in open(ioe):
            if line.startswith("seqname"):
                continue
            chrom, gene, eidx, T1, T2, = line.strip(
                "\n"
            ).split("\t")
            if chrom not in events:
                events[chrom] = {}
            if gene not in events[chrom]:
                events[chrom][gene] = []
            T1 = set(T1.split(","))
            T2 = set(T2.split(","))
            assert len(T1) < len(T2)
            T2 = T2 - T1
            events[chrom][gene].append((eidx, T1, T2))
    return events


def select_events(events, P):
    selected_events = {}
    for chrom in events:
        selected_events[chrom] = {}
        for gene in events[chrom]:
            if random.randint(1, 100) <= P:
                event = random.choice(events[chrom][gene])
                selected_events[chrom][gene] = event
    return selected_events


def split_annotation(gtf, WD, events):
    events_log = open(os.path.join(WD, "events.txt"), "w")
    anngtf_path = os.path.join(WD, "annotation.gtf")
    simgtf_path = os.path.join(WD, "simulation.gtf")
    anngtf = open(anngtf_path, "w")
    simgtf = open(simgtf_path, "w")
    for gene in gtf.features_of_type("gene"):
        anngtf.write(str(gene) + "\n")
        simgtf.write(str(gene) + "\n")
        chrom = gene.seqid
        gidx = gene.id
        annotated_transcripts = set()
        simulation_transcripts = set()
        if gidx in events[chrom]:
            event = events[chrom][gidx]
            eidx = event[0]
            _etype, chrom, *positions, strand = eidx.split(";")[1].split(":")

            T1, T2 = "", ""
            if _etype == "SE":
                ab, cd = positions
                intron1 = [int(x) for x in ab.split("-")]
                intron2 = [int(x) for x in cd.split("-")]
                T1 = f"{intron1[0] + 1}-{intron1[1] - 1}|{intron1[1]}-{intron2[0]}|{intron2[0] + 1}-{intron2[1] - 1}"
                T2 = f"{intron1[0] + 1}-{intron2[1] - 1}"
            elif _etype == "A3" or _etype == "A5":
                ab, cd = positions
                intron1 = [int(x) for x in ab.split("-")]
                intron2 = [int(x) for x in cd.split("-")]
                # we want shorter intron first
                if abs(intron1[1] - intron1[0]) > abs(intron1[1] - intron1[0]):
                    intron1, intron2 = intron2, intron1
                T1 = f"{intron1[0] + 1}-{intron1[1] - 1}"
                T2 = f"{intron2[0] + 1}-{intron2[1] - 1}"
            elif _etype == "RI":
                a, bc, d = positions
                a, d = int(a), int(d)
                intron1 = [int(x) for x in bc.split("-")]
                T1 = f"{a}-{d}"
                T2 = (
                    f"{a}-{intron1[0]}|{intron1[0]+1}-{intron1[1] - 1}|{intron1[1]}-{d}"
                )
            elif _etype == "MX":
                ab, cd, ef, gh = positions
                intron1 = [int(x) for x in ab.split("-")]
                intron2 = [int(x) for x in cd.split("-")]
                intron3 = [int(x) for x in ef.split("-")]
                intron4 = [int(x) for x in gh.split("-")]
                T1 = f"{intron1[0]+1}-{intron1[1]-1}|{intron1[1]}-{intron2[0]}|{intron2[0]+1}-{intron2[1]-1}"
                T2 = f"{intron3[0]+1}-{intron3[1]-1}|{intron3[1]}-{intron4[0]}|{intron4[0]+1}-{intron4[1]-1}"

            annotated_transcripts = event[1]
            simulation_transcripts = event[2]
            i = random.randint(1, 2)
            if i == 2:
                annotated_transcripts, simulation_transcripts = (
                    simulation_transcripts,
                    annotated_transcripts,
                )
                T2, T1 = T1, T2
                if _etype == "SE":
                    _etype = "CE"
                elif _etype in ["A3", "A5"]:
                    _etype += "'"
                elif _etype == "RI":
                    _etype += "'"

            print(
                chrom,
                gidx,
                _etype,
                event[0],
                T1,
                T2,
                ",".join(annotated_transcripts),
                ",".join(simulation_transcripts),
                sep="\t",
                file=events_log,
            )

        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            tidx = transcript.id
            anno_flag = len(annotated_transcripts) == 0 or tidx in annotated_transcripts
            sim_flag = (
                len(simulation_transcripts) == 0 or tidx in simulation_transcripts
            )
            # if both are false, then transcript is not part of the event
            if not anno_flag and not sim_flag:
                anno_flag = True
                sim_flag = True
            if anno_flag:
                anngtf.write(str(transcript) + "\n")
            if sim_flag:
                simgtf.write(str(transcript) + "\n")
            for feature in gtf.children(transcript, order_by="start"):
                if anno_flag:
                    anngtf.write(str(feature) + "\n")
                if sim_flag:
                    simgtf.write(str(feature) + "\n")
    anngtf.close()
    simgtf.close()
    events_log.close()

    return simgtf_path, anngtf_path
