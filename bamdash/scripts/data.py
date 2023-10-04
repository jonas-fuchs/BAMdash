"""
contains defs for data analysis
"""
# BUILT INS
import statistics

# LIBS
import pandas as pd
from Bio import SeqIO
import pysam
from pysam import VariantFile


def make_stat_substring(stat_string, name, value):
    return stat_string + f"<b>{name}:</b> {value}    "


def make_title_string(parsed_bam, coverage_df, reference):
    """
    :param parsed_bam: parsed bam
    :param reference:reference id
    :param coverage_df: df with computed coverage
    :return: string for header
    """
    # get bam stats for correct chrom
    bam_stats = parsed_bam.get_index_statistics()[0]
    # format title string
    stat_string = ""
    stat_string = make_stat_substring(stat_string, "reference", bam_stats[0])
    stat_string = make_stat_substring(stat_string, "reference length", f"{parsed_bam.get_reference_length(reference)} bp")
    gc_content = round((sum(coverage_df["C"])+sum(coverage_df["G"]))/len(coverage_df), 2)
    for bam_stat, stat_type in zip(bam_stats[1:], ["mapped", "unmapped", "total reads"]):
        stat_string = make_stat_substring(stat_string, stat_type, bam_stat)
    stat_string = make_stat_substring(stat_string, "gc content", f"{gc_content}%")

    return stat_string


def bam_to_coverage_df(bam_file, ref):
    """
    :param bam_file: bam location
    :param ref: chrom identifier
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
    return coverage_df, make_title_string(bam, coverage_df, ref)


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
    variant_dict = {x: [] for x in ["position", "reference", "mutation", "type"] + vcf_info}
    # create vcf dictionary
    for rec in vcf.fetch():
        if rec.chrom != ref:
            continue
        variant_dict["position"].append(rec.pos)
        variant_dict["reference"].append(rec.ref)
        variant_dict["mutation"].append(rec.alts[0])
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
        if all([x != None for x in variant_dict[key]]):
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


def genbank_to_dict(gb_file, coverage_df, ref):
    """
    parses genbank to dic and computes coverage for each annotation
    :param gb_file: genbank record location
    :param coverage_df: df with computed coverages
    :param ref: chrom id
    :return: feature_dict: dictionary with all features
    """

    feature_dict = {}

    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        if gb_record.id != ref and gb_record.name != ref:
            return feature_dict
        for feature in gb_record.features:
            if feature.type not in feature_dict:
                feature_dict[feature.type] = {}
            start, stop = feature.location.start + 1, feature.location.end
            feature_dict[feature.type][f"{start} {stop}"] = {}
            feature_dict[feature.type][f"{start} {stop}"]["mean coverage"] = round(
                statistics.mean(
                    coverage_df[(coverage_df["position"] >= start) & (coverage_df["position"] <= stop)]["coverage"]
                )
            )
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

    return define_track_position(feature_dict)


def bed_to_dict(bed_file, coverage_df, ref):
    """
    parses bed file to dic and computes coverage for each annotation
    :param bed_file: bed file location
    :param coverage_df: df with computed coverages
    :param ref: chrom id
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
        # check for additional info
        if len(line) > 3:
            for element, classifier in zip(line[3:], possible_classifiers):
                bed_dict["bed annotations"][f"{start} {stop}"][classifier] = element
        # compute mean coverage
        bed_dict["bed annotations"][f"{start} {stop}"]["mean coverage"] = round(
            statistics.mean(
                coverage_df[(coverage_df["position"] >= start) & (coverage_df["position"] <= stop)]["coverage"]
            )
        )
    return define_track_position(bed_dict)

