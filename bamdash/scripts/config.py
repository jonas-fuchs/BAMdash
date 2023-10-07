"""
contains config setting for BAMDASH to customize the plot
"""

# pdf settings
show_log = True

# overall layout
vcf_track_proportion = 0.3
gb_track_proportion = 0.5
bed_track_proportion = 0.2
plot_spacing = 0.05
font = "Arial"
font_size = 16

# coverage customize
coverage_fill_color = "rgba(255, 212, 135, 0.4)"
coverage_line_color = "rgba(224, 168, 68, 1)"
average_line_color = "grey"
average_line_width = 1

# track customize
track_color_scheme = "agsunset"  # for mutiple annotations tracks (genebank)
track_color_single = "rgb(145, 145, 145)"  # for single tracks (any rgb value, but no named colors)
strand_types = ["triangle-right", "triangle-left", "diamond-wide"]  # +, -, undefined strand
strand_marker_size = 8
strand_marker_line_width = 1
strand_marker_line_color = "rgba(0, 0, 0, 0.2)"
box_bed_alpha = [0.7, 0.7]  # alpha values for boxes (bed)
box_bed_size = [0.4, 0.4]  # size values for boxes (bed)
box_gb_alpha = [0.7, 0.8]  # alpha values for boxes (gb)
box_gb_size = [0.4, 0.3]  # size values for boxes (gb)

# variant customize
variant_marker_size = 13
variant_marker_line_width = 1
variant_line_color = "black"
stem_color = "grey"
stem_width = 1
snp_color = "grey"
ins_color = "blue"
del_color = "red"
