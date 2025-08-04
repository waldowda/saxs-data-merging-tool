# saxs-data-merging-tool
Tool to merge SAXS data from Xenocs X-Ray Scattering Instrument.

A PyQt5 GUI application for merging Small Angle X-ray Scattering (SAXS) datasets with overlapping Q ranges. 
This tool was mainly developed using Claude.AI. As allows, make sure the data merging is making sense for
your work.

## Features

- Interactive file selection
- Automatic merge region optimization
- Optional intensity scaling
- Real-time preview of merge regions
- Comprehensive data export with metadata
- Multi-threaded processing

## Installation

1. Clone this repository:
```bash
git clone https://github.com/YOUR_USERNAME/saxs-data-merging-tool.git
cd saxs-data-merging-tool

2. Instdall required packages
pip install -r requirements.txt

3. Run the application
You may need to set the python file to execute.
chmod +x saxs_merge_giu.py
python saxs_merge_gui

## Using the program
- enter your lower  q data in dataset 1.
- Enter your higher q data in dataset 2.
- Try the automated merge choice and evaluate it.
- Manuallly enter a q range if you prefer a different range.
- Click Merge and then Save and choose a filename / location.
- Click new analysis if you want to process another pair.

### Citation
Waldow, D. (2025). SAXS Data Merging Tool (Version 1.0.0) [Computer software]. 	Zenodo. https://doi.org/10.5281/zenodo.1673402.
