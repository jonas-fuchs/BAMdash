"""
contains defs for data analysis
"""

# LIBS
import pandas as pd


def bam_to_coverage_df(bam, ref):
    """
    :param bam: pysam bam file
    :param ref: chrom identifier
    :return: dataframe containing coverage and percentage nucleotides per position
    """
    coverage, position = [], []
    # count coverage at each pos
    coverage_base = bam.count_coverage(ref)
    # extratc overall coverage and pos
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

