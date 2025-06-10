from setuptools import setup, find_packages
from bamdash import __version__, _program

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='bamdash',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    python_requires=">=3.9",
    license_files=('LICENSE'),
    packages=find_packages(),
    install_requires=[
        "kaleido>=0.2.1",
        "pandas>=1.4.4",
        "plotly>=5.17.0",
        "pysam>=0.21.0",
        "biopython>=1.79"
    ],
    description='create a interactive coverage plot dashboard from bam files and add gb, vcf and bed tracks',
    url='https://github.com/jonas-fuchs/BAMdash',
    author='Dr. Jonas Fuchs',
    author_email='jonas.fuchs@uniklinik-freiburg.de',
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)"
    ],
    entry_points="""
    [console_scripts]
    {program} = bamdash.command:main
    """.format(program=_program),
    include_package_data=True,
    keywords=[],
    zip_safe=False
)