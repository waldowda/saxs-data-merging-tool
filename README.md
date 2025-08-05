# SAXS Data Merging Tool

A PyQt5 GUI application for merging Small Angle X-ray Scattering (SAXS) datasets with overlapping Q ranges, specifically developed for data from Xenocs X-Ray Scattering Instruments.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16734022.svg)](https://doi.org/10.5281/zenodo.16734022)

> **Note:** This tool was developed using Claude.AI. Always verify that the data merging results make sense for your specific research application.

## Features

- **Interactive file selection** with automatic format detection
- **Automatic merge region optimization** using RMS minimization
- **3-column data support** for Q, Intensity, and Error data
- **Error-weighted scaling** to match dataset intensities using inverse variance weighting
- **Proper uncertainty propagation** in overlap regions with weighted averaging
- **Optional error bar plotting** with performance-optimized display
- **Real-time preview** of merge regions and data quality assessment
- **Comprehensive data export** with detailed metadata headers including error data
- **Multi-threaded processing** to prevent GUI freezing during analysis
- **Smart Q-range validation** with automatic dataset swapping suggestions
- **Xenocs header detection** and automatic skipping

## Installation

### Prerequisites
- Python 3.7 or higher
- Required Python packages (see requirements.txt)

### Setup Instructions

1. **Clone this repository:**
   ```bash
   git clone https://github.com/waldowda/saxs-data-merging-tool.git
   cd saxs-data-merging-tool
   ```

2. **Install required packages:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Make the script executable (macOS/Linux):**
   ```bash
   chmod +x saxs_merge_gui.py
   ```

4. **Run the application:**
   ```bash
   python saxs_merge_gui.py
   ```

### Dependencies
- PyQt5 >= 5.12
- matplotlib >= 3.3
- numpy >= 1.19
- pandas >= 1.2
- scipy >= 1.6

## Usage Guide

### Step-by-Step Instructions

1. **Launch the application:**
   ```bash
   python saxs_merge_gui.py
   ```

2. **Load your datasets:**
   - **Dataset 1 (lower Q data):** Select your low-angle/long-distance measurement file
   - **Dataset 2 (higher Q data):** Select your wide-angle/short-distance measurement file
   
   The tool will automatically validate Q-range ordering and suggest swapping if needed.

3. **Configure merge settings:**
   - Choose whether to apply intensity scaling to match datasets
   - Use the "Find Auto Merge Region" button for automatic optimization
   - Or manually adjust the Q start/end values for custom merge ranges

4. **Preview and evaluate:**
   - Click "Preview Merge Region" to visualize the overlap
   - Assess data quality and noise levels in the proposed merge region
   - Adjust merge boundaries if necessary

5. **Merge datasets:**
   - Click "Merge Datasets" to perform the combination
   - Review the final merged plot showing original and combined data

6. **Save results:**
   - Click "Save Results" and choose filename/location
   - Both data file (.dat) and plot (.png) are automatically saved
   - Comprehensive metadata is included in the output file header

7. **Process additional datasets:**
   - Click "New Analysis" to reset and process another pair of datasets

## Advanced Features

### Error Data Support (v1.0.2+)
- **Automatic detection** of 3-column data files (Q, I, σ)
- **Error-weighted merging** using inverse variance weighting
- **Uncertainty propagation** throughout the merge process
- **Error bar visualization** with optional display control
- **Publication-ready output** with proper error representation

### Xenocs Instrument Support
- **Automatic header detection** and skipping for Xenocs files
- **GUI logging** of header processing status
- **Robust data parsing** with multiple fallback methods

## File Formats

Supports common SAXS data formats:
- `.dat`, `.txt` (space/tab delimited) 
- `.csv` (comma separated)
- `.xy`, `.chi` (two-column formats)

**Expected format:** 
- **2-column:** Q values, Intensity
- **3-column:** Q values, Intensity, Error (recommended)
- Automatic detection of data columns and Xenocs headers
- Error estimation from Poisson statistics when error column unavailable

## Algorithm Details

The tool uses advanced uncertainty propagation and optimization techniques:

### Merge Region Detection
- **RMS minimization** to find optimal merge regions
- **Edge avoidance** - searches central 60% of overlap to avoid noisy boundaries
- **Median-based scaling** factors for robustness against outliers
- **Common Q-grid interpolation** for accurate comparison

### Error Handling and Propagation
- **3-column support** with automatic error detection
- **Inverse variance weighting** for error-weighted scaling: `w = 1/σ²`
- **Proper uncertainty propagation** in overlap regions: `σ_merged = √(1/(w₁ + w₂))`
- **Poisson estimation** when error data unavailable: `σ = √I`

### Performance Optimizations
- **Multi-threaded processing** for merge operations
- **Error bar thinning** for large datasets in plots
- **Efficient interpolation** algorithms for smooth merging

## Author

**Dean Waldow**  
Pacific Lutheran University  
Email: waldowda@plu.edu

## Citation

If you use this tool in your research, please cite:

> Waldow, D. (2025). SAXS Data Merging Tool (Version 1.0.2) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.16734022

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Changelog

### Version 1.0.2 (Latest)
- ✅ **3-column data support** (Q, I, σ) with automatic detection
- ✅ **Error-weighted scaling** using inverse variance weighting  
- ✅ **Proper uncertainty propagation** in merge overlap regions
- ✅ **Optional error bar plotting** with GUI checkbox control
- ✅ **Enhanced user feedback** for error data detection
- ✅ **Updated save format** to include error column
- ✅ **Performance optimized** error bar display for large datasets

### Version 1.0.1
- 🔧 **Fixed syntax errors** in data loading functions
- 🔧 **Full Q-range plotting** in preview windows (no more cropping)
- 🔧 **Xenocs header GUI logging** for better user feedback
- 🔧 **Improved error handling** and user experience

### Version 1.0.0
- 🎉 **Initial release** with PyQt5 GUI interface
- 🎉 **Automatic merge region detection** using RMS optimization
- 🎉 **Optional intensity scaling** functionality
- 🎉 **Interactive plotting** and preview capabilities

## Troubleshooting

### Common Issues

**PyQt5 installation problems:**
```bash
pip install --user PyQt5  # User installation
# or
conda install pyqt        # Using conda
```

**File permission errors:**
- Ensure Python has access to your data files
- On macOS: Check System Preferences > Security & Privacy > Files and Folders

**Memory issues with large datasets:**
- Consider downsampling very large files before merging
- Close unnecessary applications to free memory

## Contributing

Bug reports, feature requests, and contributions are welcome! Please open an issue or submit a pull request on GitHub.

---

**Disclaimer:** Always validate merged results against your experimental expectations and known standards.