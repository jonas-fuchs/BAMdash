"""
contains config setting for BAMcov
"""
# gives the track proportion relative to the coverage plot
vcf_track_proportion = 0.3
gb_track_proportion = 0.5
bed_track_proportion = 0.3

# colors
coverage_fill_color = "rgba(255, 212, 135, 0.2)"  # coverage plot
coverage_line_color = "rgba(224, 168, 68, 1)"  # coverage plot
average_line_color = "grey"  # coverage plot
snp_color = "grey"  # vcf plot
ins_color = "blue"  # vcf plot
del_color = "red"  # vcf plot
stem_color = "grey"  # vcf plot
color_scheme = "agsunset"  # gb plot
box_alpha = [0.6, 0.8]  # alpha values for boxes

# size
box_size = [0.4, 0.3]  # gb plot

# markers
strand_types = ["arrow-bar-right", "arrow-bar-left", "diamond-wide"]  # +, -, undefined strand
