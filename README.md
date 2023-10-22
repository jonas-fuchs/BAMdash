
<img src="./bamdash.png" alt="bamdash" />

[![language](https://img.shields.io/badge/python-%3E3.9-green)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/github/license/jonas-fuchs/bamdash)](https://www.gnu.org/licenses/gpl-3.0)
![Static Badge](https://img.shields.io/badge/platform-linux_osx-blue)
[![DOI](https://zenodo.org/badge/700952196.svg)](https://zenodo.org/badge/latestdoi/700952196)
[![PiPy](https://img.shields.io/pypi/v/bamdash?label=pypi%20version)](https://pypi.org/project/bamdash/)
[![Downloads](https://static.pepy.tech/badge/bamdash)](https://pypi.org/project/bamdash/)

## Overview

**BAMdash lets you create interactive coverage plots from your bam file with [`plotly`](https://plotly.com/)**

- **requires** only a `.bam`, `.bai` and the reference id to which the reads where mapped
- **create** a interactive `html` for data exploration
- **create** a static image (`jpg`, `png`, `pdf`, `svg`) ready for publication
- **add** additional tracks (supported: `.vcf`, `.gb`, `.bed`)
- **annotate** tracks with additional information
- **export** annoated track data as tabular files (`.bed`, `.vcf`) or json (`.gb`)
- **developed** for viral genomics
- **customize** all plotting parameters

**Feel free to report any bugs or request new features as issues!**


## Automatic annotation

BAMdash automatically computes serveral statistics:

- if `-bs` is > 1 it computes the mean over the bin size in the coverage plot
- for each track it computes recovery and mean coverage (set `-c` for the min coverage) for each element in the track
- if a `*.vcf` is provided it annotates `TRANSITION`/`TRANSVERSION` and type of exchange (`SNP`, `DEL`, `INS)

If a `*.gb`and `*.vcf` is provided BAMdash computes the aminoacid exchange and the effect in the CDS (inspired by but not as powerful as [snpeff](http://pcingola.github.io/SnpEff/snpeff)). SNP and INDEL vcf annotation supports:

- `START_LOST`: INDEL or SNP start at the CDS and result in a start loss
- `STOP_LOST`: INDEL or SNP result in the loss of the stop codon
- `STOP_GAINED`: INDEL or SNP result in an additional stop codon
- `SYN`: SNP does not lead to an amino acid change
- `NON-SYN`: SNP leads to an amino acid change 
- `AC_INSERTION`: INS that does not change already present amino acids
- `AC_CHANGE+AC_INSERTION`: INS where the affected codon is also non-syn
- `AC_DELETION`: DEL that does not change already present amino acids
- `AC_CHANGE+AC_DELETION`: DEL where the affected codon is also non-syn

The nomenclature for the aminoacid effect is pretty simplified:

- `A58Y` - Exchange at pos 58 from A to Y
- `A58YY`- Exchange at pos 58 from A to Y and insertion of an additional Y
- `FA58Y`- Exchange at pos 58 from A to Y and deletion of the prior F
- `A58fsX` - Frameshift at pos 58

## Example
<img src="./example.gif" alt="example" />

## Installation

### via pip (recommened):
```shell
pip install bamdash
```
### from this repo:
```shell
git clone https://github.com/jonas-fuchs/BAMdash
cd BAMcov
```
and then install BAMdash with:
```shell
pip install -r requirements.txt
```
or:
```shell
pip install .
```
That was already it. To check if it worked:

```shell
bamdash -v
```
You should see the current BAMdash version.

## Usage

```shell
usage: 	

bamdash -b bam_file_path -r reference_id [additional arguments]
```
```
full usage:

  -h, --help            show this help message and exit
  -b  , --bam           bam file location
  -r  , --reference     seq reference id
  -bs  , --binsize      bins for the coverage plot
  -t [track_1 ...], --tracks [track_1 ...]
                        file location of tracks
  -c 5, --coverage 5    minimum coverage
  --slider, --no-slider
                        show slider (default: False)
  -e None, --export_static None
                        export as png, jpg, pdf, svg
  -d px px, --dimensions px px
                        width and height of the static image in px
  --dump, --no-dump     dump annotated track data (default: False)
  -v, --version         show program's version number and exit
```

## Cutomization

BAMcov plotting settings can be adjusted in in the `config.py`. Therefore, you have to clone this repo.

Go to the configs location:
```shell
cd BAMdash/bamdash/scripts/
```
And open the `config.py` with a text editor, e.g.:
```shell
gedit config.py
```
and adjust the settings:
```python
# pdf settings
show_log = True

# overall layout
vcf_track_proportion = 0.3
gb_track_proportion = 0.5
bed_track_proportion = 0.2
plot_spacing = 0.05

# coverage customize
coverage_fill_color = "rgba(255, 212, 135, 0.2)"
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
box_bed_alpha = [0.6, 0.6]  # alpha values for boxes (bed)
box_bed_size = [0.4, 0.4]  # size values for boxes (bed)
box_gb_alpha = [0.6, 0.8]  # alpha values for boxes (gb)
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
```
To apply these new settings just repeat the installation procedure in the BAMdash dir:
```shell
pip install .
```

<a href="https://www.buymeacoffee.com/jofox" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>


---

**Important disclaimer:**
*The code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
