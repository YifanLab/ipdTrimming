<ipdTrimming: a pipeline for 6mA calling>
    Copyright (C) <2025>  <Wentao Yang>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.




# 6mA_footprint
## Here is the comparison of pipelines between ipdSummary and ipdTrimming.

![image](https://github.com/user-attachments/assets/610a723a-c0bf-4203-a37e-1d0d260944c3)

For ipdTrimming, a Local trimming of IPD outliers at the subread level was introduced into IPD conversion: for each site in a DNA molecule, the top 10% of subread IPD values were removed; the remaining subread IPD values were then averaged to generate the CCS IPD value.  This was implemented by a modified CCS module. 

```console
usc@hpc:~$ ccs -j 20 --hifi-kinetics  subread.bam hifi.bam
    
```

For this modified CCS module, we get the CCS bam file with trimmed IPD value. Following a movie-time normalization was applied to compensate for variations in the polymerase elongation rate on individual DNA molecules, each CCS IPD value was normalized against the CCS IPD value averaged across all sites in a DNA molecule.  This CCS IPD value was then compared to that of an unmodified base with the same sequence context, provided by a kmer model pretrained on the Sequel data. The last two steps were implemented by the script ipdRatiocalculator_FromCCS.py 

```console
usc@hpc:~$ python ipdRatiocalculator_FromCCS.py hifi.bam hifi.withIPDr.bam
#the cpu number can be set in this script.
```

Extract the IPDr for all A sites in individual molecule
```console
usc@hpc:~$ python bamextractallAx_IPDvalue.py hifi.withIPDr.bam hifi.withIPDr_allA.xls
#the effective coverage (ec) can be set in this script, the default value is 20.
```

## Tutorial for ipdRatiocalculator_fromCCS.py
ipdRatiocalculator_fromCCS.py is a Python script for processing HiFi kinetic BAM files to compute IPD (inter-pulse duration) ratios using a pre-trained model.

### Prerequisites

Before running the script, ensure that the following dependencies are installed in a Python 3 environment:

### 1.1 Install Dependencies via Conda
It is recommended to use Conda to manage dependencies:
```console
conda create -n ipdcalc_env python=3.7 -y
conda activate ipdcalc_env
conda install -c bioconda numpy pandas tqdm pysam pbcore
```
### 1.2 Install kineticsTools
```console
git clone https://github.com/PacificBiosciences/kineticsTools.git
cd kineticsTools
pip install .
```

### 1.3 Required External Software
PacBio SMRTLink Suite (Provides the kinetics tools)

### Preparing Input Files
The script requires:

    A BAM file: Unaligned HiFi kinetic BAM file processed with ccs-kinetics-bystrandify
    A lookup table (SP3-C3.npz.gz): This file is hardcoded into the script. Ensure the correct path is used. If you don't find this file, you can download from this site: https://github.com/PacificBiosciences/kineticsTools/blob/master/kineticsTools/resources/SP3-C3.npz.gz.

### Troubleshooting
Issue: Missing Dependencies

If you encounter an error like:
```console
ModuleNotFoundError: No module named 'pbcore'
```
Ensure dependencies are installed using Conda or Pip:
```console
conda install -c bioconda pbcore
pip install pbcore
```


