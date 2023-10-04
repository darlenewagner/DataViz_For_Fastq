# DataViz_For_Fastq
A showcase of functions using Biopython, FastqGeneralIterator, pandas, and matplotlib for teaching Big Data parsing for NGS

## INSTALLATION

#### 1. Getting started after GIT download
`cd DataViz_For_Fastq`

`gunzip Life_Expectancy_US/US_A_USALEEP.csv.gz` 

#### 2. Install Python 3.9.X or higher
Alternatively, in a multi-user environment, load Python module:

`module load Python/3.9.1`

#### 3. Set up Python virtual environment in current folder
Exact path to bin/python may differ according to Python version.

`virtualenv -p /apps/x86_64/python/3.9.1/bin/python ./`

#### 4. Install prerequisite Python packages
`bin/pip install biopython`

`bin/pip install matplotlib`

`bin/pip install pandas`

`bin/pip install mpld3`

`bin/pip install plotly`

## USAGE EXAMPLES

#### 1. Plotting NGS reads lengths from four .fastq files
`bin/python BioPython_LengthViolinFastq_Reads.py Viral_MiSeq/SC2_X.R1.fastq Viral_MiSeq/SC2_X.R2.fastq Viral_iSeq/SC2_Y.R1.fastq Viral_iSeq/SC2_Y.R2.fastq --outputType P`

#### 2. Plotting NGS reads PHRED (quality) scores from two fastq files
`BioPython_PHREDViolinFastq_Reads.py Viral_MiSeq/SC2_X.R1.fastq Viral_MiSeq/SC2_X.R2.fastq --outputType P --titleString "Raw Reads PHRED Scores" --showQ30 Y`

Note: For faster runtime, `--showQ30 N`

#### 3. Plotting Census Tract Life Expectancies for States from USALEEP Dataset
`bin/python Pandas_ParseForViolin_CSV.py US_A_USALEEP.csv --outputType P --stateCodes 6 56 --titleString "U.S., California, and Wyoming Life Expectancy,"`

#### 4. Interactive Plots of Life Expectancies for States from USALEEP Dataset
`bin/python Plotly_ParseForViolin_to_HTML.py Life_Expectancy_US/US_A_USALEEP.csv --stateCodes 15 6 30 --titleString "U.S., Hawaii, California, and Montana"`

Note: For Plotly_ParseForViolin_to_HTML.py, output path should be set to a folder that enables point-and-click opening of the output HTML.
