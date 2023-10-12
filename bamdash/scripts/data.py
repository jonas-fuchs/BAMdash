"""
contains defs for data analysis
"""
# BUILT INS
import statistics

# LIBS
import pandas as pd
from Bio import Seq
from Bio import SeqIO
import pysam
from pysam import VariantFile


def make_stat_substring(stat_string, name, value):
    return stat_string + f"<b>{name}:</b> {value}    "


def get_coverage_stats(coverage_df, start, stop, min_cov):
    """
    :param coverage_df:
    :param start: start of an annotation
    :param stop: stop of an annotation
    :param min_cov: min coverage to consider covered
    :return: stats
    """
    df_subset = coverage_df[(coverage_df["position"] >= start) &
                    (coverage_df["position"] <= stop) &
                    (coverage_df["coverage"] > min_cov)]
    if df_subset.empty:
        mean_coverage = 0
        recovery = 0
    else:
        mean_coverage = statistics.mean(df_subset["coverage"])
        recovery = len(df_subset["coverage"])/(stop-start+1)*100

    return round(mean_coverage), round(recovery, 2)


def make_title_string(parsed_bam, coverage_df, reference, min_cov):
    """
    :param parsed_bam: parsed bam
    :param reference:reference id
    :param coverage_df: df with computed coverage
    :param min_cov: min coverage to consider covered
    :return: string for header
    """
    # get bam stats for correct chrom
    bam_stats = parsed_bam.get_index_statistics()[0]
    mean, rec = get_coverage_stats(coverage_df, min(coverage_df["position"]), max(coverage_df["position"]), min_cov)
    # format title string
    stat_string = ""
    stat_string = make_stat_substring(stat_string, "reference", bam_stats[0])
    stat_string = make_stat_substring(stat_string, "reference length", f"{parsed_bam.get_reference_length(reference)} bp")
    gc_content = round((sum(coverage_df["C"])+sum(coverage_df["G"]))/len(coverage_df), 2)
    for bam_stat, stat_type in zip(bam_stats[1:], ["mapped", "unmapped", "total"]):
        stat_string = make_stat_substring(stat_string, stat_type, bam_stat)
    stat_string = make_stat_substring(stat_string, "<br>mean coverage", mean)
    stat_string = make_stat_substring(stat_string, "recovery", rec)
    stat_string = make_stat_substring(stat_string, "gc content", f"{gc_content}%")

    return stat_string


def bam_to_coverage_df(bam_file, ref, min_cov):
    """
    :param bam_file: bam location
    :param ref: chrom identifier
    :param min_cov: min coverage to consider covered
    :return: dataframe containing coverage and percentage nucleotides per position and title with stats
    """
    # parse bam
    bam = pysam.AlignmentFile(bam_file, "rb")

    coverage, position = [], []
    # count coverage at each pos
    coverage_base = bam.count_coverage(ref)
    # extract overall coverage and pos
    for index, [A_count, C_count, G_count, T_count] in enumerate(
            zip(coverage_base[0], coverage_base[1], coverage_base[2], coverage_base[3])):
        coverage.append(A_count + C_count + G_count + T_count)
        position.append(index + 1)
    # return the coverage dataframe
    coverage_df = pd.DataFrame(
        list(zip(
            position,
            coverage,
            [round(x / y * 100, 1) if y != 0 else 0 for x, y in zip(coverage_base[0], coverage)],
            [round(x / y * 100, 1) if y != 0 else 0 for x, y in zip(coverage_base[1], coverage)],
            [round(x / y * 100, 1) if y != 0 else 0 for x, y in zip(coverage_base[2], coverage)],
            [round(x / y * 100, 1) if y != 0 else 0 for x, y in zip(coverage_base[3], coverage)],
        )),
        columns=[
            "position",
            "coverage",
            "A",
            "C",
            "G",
            "T"
        ])
    return coverage_df, make_title_string(bam, coverage_df, ref, min_cov)


def vcf_to_df(vcf_file, ref):
    """
    :param vcf_file: vcf location
    :param ref: chrom id
    :return: vcf dictionary
    """
    # parse vcf
    vcf = VariantFile(vcf_file)

    # extract vcf info
    vcf_info = list(vcf.header.info)
    variant_dict = {x: [] for x in ["position", "reference", "mutation", "type", "point mutation type"] + vcf_info}
    # create vcf dictionary
    for rec in vcf.fetch():
        if rec.chrom != ref:
            continue
        variant_dict["position"].append(rec.pos)
        variant_dict["reference"].append(rec.ref)
        variant_dict["mutation"].append(rec.alts[0])
        # annotate the point mutation type
        base_exchange = f"{rec.ref}{rec.alts[0]}"
        if len(base_exchange) == 2:
            if base_exchange in ["AG", "GA", "CT", "TC"]:
                variant_dict["point mutation type"].append("TRANSITION")
            else:
                variant_dict["point mutation type"].append("TRANSVERSION")
        else:
            variant_dict["point mutation type"].append("NONE")

        # get mutation type
        if len(rec.alts[0]) > len(rec.ref):
            variant_dict["type"].append("INS")
        elif len(rec.alts[0]) < len(rec.ref):
            variant_dict["type"].append("DEL")
        else:
            variant_dict["type"].append("SNP")
        # populate table with info fields
        if vcf_info:
            for info_field in vcf_info:
                if info_field in rec.info:
                    variant_dict[info_field].append(rec.info[info_field])
                else:
                    variant_dict[info_field].append(None)
    # clear all keys that have None appended
    empty_keys = []
    for key in variant_dict:
        if all([x is not None for x in variant_dict[key]]):
            continue
        empty_keys.append(key)
    if empty_keys:
        [variant_dict.pop(key) for key in empty_keys]

    return pd.DataFrame.from_dict(variant_dict)


def define_track_position(feature_dict):
    """
    adds track info to feature dictionary
    :param feature_dict: feature dict without track info
    :return: feature dict with track info
    """
    # add the track (y position) to plot the feature in
    # remember for each track the largest stop
    new_track = 0
    for feature_type in feature_dict:
        track_stops = [0]
        # check if a start of a gene is smaller than the stop of the current track
        # -> move to new track
        for annotation in feature_dict[feature_type]:
            track = 0
            positions = [int(x) for x in annotation.split(" ")]
            while positions[0] < track_stops[track]:
                track += 1
                # if all prior tracks are potentially causing an overlap
                # create a new track and break
                if len(track_stops) <= track:
                    track_stops.append(0)
                    break
            # in the current track remember the stop of the current annotation
            track_stops[track] = positions[1]
            # and indicate the track in the dict
            feature_dict[feature_type][annotation]["track"] = track + new_track
            # for each annotation make sure to always add another track
        new_track += len(track_stops)

    return feature_dict


def genbank_to_dict(gb_file, coverage_df, ref, min_cov):
    """
    parses genbank to dic and computes coverage for each annotation
    :param gb_file: genbank record location
    :param coverage_df: df with computed coverages
    :param ref: chrom id
    :param min_cov: min coverage to consider covered
    :return: feature_dict: dictionary with all features
    """

    feature_dict = {}
    seq = ""

    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        if gb_record.id != ref and gb_record.name != ref:
            break
        seq = gb_record.seq
        for feature in gb_record.features:
            if feature.type not in feature_dict:
                feature_dict[feature.type] = {}
            start, stop = feature.location.start + 1, feature.location.end
            feature_dict[feature.type][f"{start} {stop}"] = {}
            mean_cov, rec = get_coverage_stats(coverage_df, start, stop, min_cov)
            feature_dict[feature.type][f"{start} {stop}"]["start"] = start
            feature_dict[feature.type][f"{start} {stop}"]["stop"] = stop
            feature_dict[feature.type][f"{start} {stop}"]["mean coverage"] = mean_cov
            feature_dict[feature.type][f"{start} {stop}"]["% recovery"] = rec
            # define strand info
            if feature.strand == 1:
                feature_dict[feature.type][f"{start} {stop}"]["strand"] = "+"
            elif feature.strand == -1:
                feature_dict[feature.type][f"{start} {stop}"]["strand"] = "-"
            else:
                feature_dict[feature.type][f"{start} {stop}"]["strand"] = "NA"
            # populate dictionary with feature infos
            for qualifier in feature.qualifiers:
                feature_dict[feature.type][f"{start} {stop}"][qualifier] = feature.qualifiers[qualifier][0]

    if feature_dict:
        return define_track_position(feature_dict), seq
    else:
        return feature_dict, seq


def bed_to_dict(bed_file, coverage_df, ref, min_cov):
    """
    parses bed file to dic and computes coverage for each annotation
    :param bed_file: bed file location
    :param coverage_df: df with computed coverages
    :param ref: chrom id
    :param min_cov: min coverage to consider covered
    :return: bed_dict: dictionary with all features
    """
    # first extract as list of list to be able to sort
    # (otherwise track info will be incorrect)
    bed_line_list = []
    for line in open(bed_file, "r"):
        if line.startswith('#'):
            continue
        line_elements = line.strip().split("\t")
        # bed file has at least chrom, start, stop
        if len(line_elements) < 3:
            continue
        # convert start, stop to integers
        line_elements[1], line_elements[2] = int(line_elements[1])+1, int(line_elements[2])
        if line_elements[0] != ref:
            continue
        bed_line_list.append(line_elements)
    # sort by start
    bed_line_list = sorted(bed_line_list, key=lambda x: x[1])

    # populate dictionary
    bed_dict = {"bed annotations": {}}
    possible_classifiers = ["name", "score", "strand"]
    for line in bed_line_list:
        start, stop = int(line[1]), int(line[2])
        bed_dict["bed annotations"][f"{start} {stop}"] = {}
        bed_dict["bed annotations"][f"{start} {stop}"]["start"] = start
        bed_dict["bed annotations"][f"{start} {stop}"]["stop"] = stop
        # check for additional info
        if len(line) > 3:
            for element, classifier in zip(line[3:], possible_classifiers):
                bed_dict["bed annotations"][f"{start} {stop}"][classifier] = element
        # compute mean coverage
        mean_cov, rec = get_coverage_stats(coverage_df, start, stop, min_cov)
        bed_dict["bed annotations"][f"{start} {stop}"]["mean coverage"] = mean_cov
        bed_dict["bed annotations"][f"{start} {stop}"]["recovery"] = rec

    return define_track_position(bed_dict)


def annotate_vcf_df(vcf_df, cds_dict, seq):
    """
    annotate mutations for their as effect
    :param vcf_df: dataframe with vcf data
    :param cds_dict: dictionary with all coding sequences
    :param seq: sequence of the reference
    :return: annotated df
    """

    # define the potential identifiers to look for
    potential_cds_identifiers = ["protein_id", "gene", "locus_tag", "product"]

    proteins, as_exchange, as_effect = [], [], []

    for variant in vcf_df.iterrows():
        proteins_temp, as_exchange_temp, as_effect_temp = [], [], []
        pos, mut_type, mut = variant[1]["position"], variant[1]["type"], variant[1]["mutation"]
        # check each cds for potential mutations
        for cds in cds_dict:
            # check if a protein identifier is present
            if not any(identifier in potential_cds_identifiers for identifier in cds_dict[cds]):
                continue
            start, stop = cds_dict[cds]["start"], cds_dict[cds]["stop"]
            # check if pos is in range
            if pos not in range(start, stop):
                continue
            # now we know the protein
            for identifier in potential_cds_identifiers:
                if identifier in cds_dict[cds]:
                    proteins_temp.append(cds_dict[cds][identifier])
                    break
            # at the moment only check SNPs
            if mut_type != "SNP":
                as_exchange_temp.append("AMBIG"), as_effect_temp.append(variant[1]["type"])
                continue
            # now we can analyse for a potential as exchange
            if "codon_start" in cds_dict[cds]:
                codon_start = int(cds_dict[cds]["codon_start"])
            else:
                codon_start = 1
            strand, seq_mut = cds_dict[cds]["strand"], str(seq)
            # mutate the sequence and get the CDS nt seq
            seq_mut = Seq.Seq(seq_mut[:pos-1] + mut + seq_mut[pos:])
            seq_cds, seq_mut_cds = seq[start+codon_start-2:stop], seq_mut[start+codon_start-2:stop]
            # translate strand depend
            if strand == "+":
                cds_trans, cds_mut_trans = seq_cds.translate(), seq_mut_cds.translate()
            else:
                cds_trans, cds_mut_trans = seq_cds.reverse_complement().translate(), seq_mut_cds.reverse_complement().translate()
            # get the mutation by searching for a diff between string  - prop a bit slow
            mut_string = [f"{x}{i+1}{y}" for i, (x, y) in enumerate(zip(cds_trans, cds_mut_trans)) if x != y]
            # now append to list
            if not mut_string:
                as_exchange_temp.append("NONE"), as_effect_temp.append("SYN")
            else:
                as_exchange_temp.append(mut_string[0]), as_effect_temp.append("NON-SYN")
        # check if we did not find a protein
        if not proteins_temp:
            proteins.append("NONE"), as_exchange.append("NONE"), as_effect.append("NONE")
        # else append all mutations found in all cds
        elif len(proteins_temp) == 1:
            proteins.append(proteins_temp[0]), as_exchange.append(as_exchange_temp[0]), as_effect.append(as_effect_temp[0])
        else:
            proteins.append(" | ".join(proteins_temp)), as_exchange.append(" | ".join(as_exchange_temp)), as_effect.append(" | ".join(as_effect_temp))

    # annotate and return df
    vcf_df["protein"], vcf_df["effect"], vcf_df["as mutation"] = proteins, as_effect, as_exchange
    return vcf_df


def annotate_vcfs_in_tracks(track_data):
    """
    annotates all vcf in the tracks
    :param track_data
    :return: track data with annotated vcfs
    """
    # annotate mutation effects in vcf df if gb is present
    index_positions = [[], []]  # index of gb and vcf
    for index, track in enumerate(track_data):
        if "gb" in track[1]:
            index_positions[0].append(index)
        elif "vcf" in track[1]:
            index_positions[1].append(index)
    # check if there is a gb file and a vcf file
    if index_positions[0] and index_positions[1]:
        # and only one gb file
        gb_indices = index_positions[0]
        if len(gb_indices) > 1:
            print("WARNING: cannot annotate from multiple *.gb files!")
            return track_data
            # annotate each vcf df
        for vcf_track in index_positions[1]:
            # check if CDS is present
            if "CDS" not in track_data[gb_indices[0]][0]:
                continue
            # check again that there is no empty data
            if not track_data[vcf_track][0].empty:
                # and then annotate the vcf dict
                track_data[vcf_track][0] = annotate_vcf_df(
                    track_data[vcf_track][0],  # vcf df
                    track_data[gb_indices[0]][0]["CDS"],  # CDS dict
                    track_data[gb_indices[0]][2]  # seq
                )

    return track_data

