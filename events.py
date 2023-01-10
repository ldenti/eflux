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
    c1gtf_path = os.path.join(WD, "condition1.gtf")
    c2gtf_path = os.path.join(WD, "condition2.gtf")
    ugtf_path = os.path.join(WD, "unused.gtf")
    c1gtf = open(c1gtf_path, "w")
    c2gtf = open(c2gtf_path, "w")
    ugtf = open(ugtf_path, "w")
    for gene in gtf.features_of_type("gene"):
        chrom = gene.seqid
        gidx = gene.id
        # Split transcripts per condition
        c1_transcripts = set()
        c2_transcripts = set()
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

            c1_transcripts = event[1]
            c2_transcripts = event[2]
            i = random.randint(1, 2)
            if i == 2:
                c1_transcripts, c2_transcripts = (
                    c2_transcripts,
                    c1_transcripts,
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
                ",".join(c1_transcripts),
                ",".join(c2_transcripts),
                sep="\t",
                file=events_log,
            )

        # Dump GTFs
        unused = len(c1_transcripts) == 0 and len(c2_transcripts) == 0
        if unused:
            # no event
            ugtf.write(str(gene) + "\n")
        else:
            c1gtf.write(str(gene) + "\n")
            c2gtf.write(str(gene) + "\n")

        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            tidx = transcript.id
            inc1, inc2 = tidx in c1_transcripts, tidx in c2_transcripts
            if unused:
                ugtf.write(str(transcript) + "\n")
            else:
                if inc1:
                    c1gtf.write(str(transcript) + "\n")
                elif inc2:
                    c2gtf.write(str(transcript) + "\n")
            for feature in gtf.children(transcript, order_by="start"):
                if unused:
                    ugtf.write(str(transcript) + "\n")
                else:
                    if inc1:
                        c1gtf.write(str(feature) + "\n")
                    elif inc2:
                        c2gtf.write(str(feature) + "\n")
    ugtf.close()
    c1gtf.close()
    c2gtf.close()
    events_log.close()

    return c1gtf_path, c2gtf_path
