"""
contains defs for data analysis
"""
# BUILT INS
import statistics

# LIBS
import pandas as pd
from Bio import SeqIO


def bam_to_coverage_df(bam, ref):
    """
    :param bam: pysam bam file
    :param ref: chrom identifier
    :return: dataframe containing coverage and percentage nucleotides per position
    """
    coverage, position = [], []
    # count coverage at each pos
    coverage_base = bam.count_coverage(ref)
    # extract overall coverage and pos
    for index, [A_count, C_count, G_count, T_count] in enumerate(
            zip(coverage_base[0], coverage_base[1], coverage_base[2], coverage_base[3])):
        coverage.append(A_count + C_count + G_count + T_count)
        position.append(index + 1)
    # return the coverage dataframe
    return pd.DataFrame(
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


def vcf_to_df(vcf, ref):
    """
    :param vcf: read vcf
    :param ref: chrom id
    :return: vcf dictionary
    """

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
    track_stops = [0]
    for feature_type in feature_dict:
        track = 0
        # check if a start of a gene is smaller than the stop of the current track
        # -> move to new track
        for annotation in feature_dict[feature_type]:
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
            feature_dict[feature_type][annotation]["track"] = track
        # for each annotation make sure to always add another track
        track += 1

    return feature_dict


def genbank_to_dict(infile, coverage_df):
    """
    parses genbank to dic and computes coverage for each annotation
    :param infile: genbank record
    :param coverage_df: df with computed coverages
    :return: feature_dict: dictionary with all features
    """

    feature_dict = {}

    for gb_record in SeqIO.parse(open(infile, "r"), "genbank"):
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
