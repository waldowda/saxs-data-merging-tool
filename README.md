# SAXS Data Merging Tool

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16734022.svg)](https://doi.org/10.5281/zenodo.16734022)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A professional graphical tool for merging Small Angle X-ray Scattering (SAXS) datasets with overlapping Q ranges. Features automatic merge region detection, optional intensity scaling, error propagation, and interactive visualization. This is written primarily through AI. As always, verify that the program functions for your needs.

## Features

### Core Functionality
- **Intelligent Data Loading**: Automatic detection and handling of various SAXS data formats including Xenocs headers
- **Automatic Q-Range Validation**: Detects when datasets are loaded in wrong order with automatic swap option
- **Smart Merge Region Detection**: Automatically finds optimal merge regions based on data overlap quality
- **Error Propagation**: Full support for uncertainty data with proper error-weighted averaging
- **Interactive Visualization**: Real-time preview of merge regions and final results

### Advanced Capabilities
- **Multi-column Support**: Handles 2-column (Q, I) or 3-column (Q, I, σ) data formats
- **Intensity Scaling**: Optional scaling with error-weighted calculation for accurate dataset matching
- **Quality Assessment**: Visual feedback for merge quality with comprehensive plotting
- **Comprehensive Export**: Saves merged data with full metadata and publication-ready plots

## Installation

### Requirements
- Python 3.7 or higher
- PyQt5 >= 5.12
- matplotlib >= 3.3
- numpy >= 1.19
- pandas >= 1.2
- scipy >= 1.6

### Install Dependencies

```bash
# Using pip
pip install PyQt5 matplotlib numpy pandas scipy

# Using conda
conda install pyqt matplotlib numpy pandas scipy
```

### Download and Run

```bash
# Clone or download the repository
git clone https://github.com/waldowda/saxs-data-merging-tool.git
cd saxs-data-merging-tool

# Run the application
python saxs_merge_gui.py
```

## Usage

### Quick Start
1. **Launch the application**: Run `python saxs_merge_gui.py`
2. **Load datasets**: 
   - Select Dataset 1 (lower Q data) using the first Browse button
   - Select Dataset 2 (higher Q data) using the second Browse button
3. **Automatic validation**: The tool will check Q-range ordering and offer to swap if needed
4. **Find merge region**: Click "Find Auto Merge Region" for optimal automatic detection
5. **Preview**: Click "Preview Merge Region" to visualize the proposed merge
6. **Merge**: Click "Merge Datasets" to perform the actual merging
7. **Save**: Use "Save Results" to export the merged data and plots

### Supported File Formats
- `.dat` - Standard SAXS data files
- `.txt` - Text files with columnar data
- `.csv` - Comma-separated values
- `.xy` - XY data format
- `.chi` - Chi format files
- **Xenocs format** - Automatically detected and handled

### Data Column Formats
- **2-column**: Q, Intensity (errors estimated from Poisson statistics)
- **3-column**: Q, Intensity, Error/Uncertainty

## Workflow

### 1. Data Loading and Validation
The tool automatically:
- Detects file headers (including Xenocs format)
- Validates Q-range ordering between datasets
- Estimates errors if not provided
- Checks for data overlap

### 2. Merge Region Optimization
- **Automatic detection**: Uses mathematical optimization to find the best merge point
- **Manual adjustment**: Fine-tune merge start/end points using spin boxes
- **Visual preview**: Real-time plotting shows merge quality

### 3. Data Merging Process
- **Scaling calculation**: Error-weighted scaling factors for accurate intensity matching
- **Overlap averaging**: Proper uncertainty propagation in merge regions
- **Quality control**: Comprehensive validation of merge results

### 4. Results and Export
- **Visualization**: Dual-panel plots showing original and merged data
- **Metadata preservation**: Complete processing history in output files
- **Publication-ready**: High-resolution plots suitable for manuscripts

## Example Output

The merged dataset includes:
```
# Merged SAXS data from dataset1.dat and dataset2.dat
# Merge range: Q = 0.085000 to 0.150000
# Scaling applied: Yes
# Scaling factor (dataset 2): 1.2345
# Error propagation: Weighted averaging in overlap region
# Total points: 1247
# Q range: 0.005000 to 2.500000
# Q(A^-1)  Intensity  Error
5.000000e-03  1.234567e+02  5.678901e+00
...
```

## Advanced Features

### Error Handling
- **Automatic error estimation**: Uses Poisson statistics when error column unavailable
- **Zero/negative error handling**: Automatically corrects problematic error values
- **Weighted calculations**: All scaling and averaging operations use proper error weighting

### Quality Control
- **Visual validation**: Multiple plot views for assessing merge quality
- **Statistical feedback**: Quantitative measures of merge success
- **Diagnostic information**: Comprehensive logging of all processing steps

### Customization Options
- **Scaling toggle**: Enable/disable intensity scaling
- **Error bar display**: Toggle error bars in plots for clarity
- **Manual merge regions**: Override automatic detection when needed

## Troubleshooting

### Common Issues

**File loading errors:**
- Ensure data files have numeric columns
- Check for consistent delimiters (space, tab, comma)
- Verify Q values are in ascending order

**No overlap detected:**
- Check that datasets actually overlap in Q-range
- Verify datasets are in correct order (Dataset 1 = lower Q)
- Consider if data units are consistent

**Poor merge quality:**
- Try adjusting merge region boundaries manually
- Check if scaling should be enabled/disabled
- Verify data quality in overlap region

**Performance issues:**
- For very large datasets (>10,000 points), consider data reduction
- Close other applications to free memory
- Use SSD storage for faster file access

## Technical Details

### Algorithms
- **Merge point optimization**: Minimizes RMS difference in overlap regions
- **Scaling calculation**: Error-weighted least squares fitting
- **Interpolation**: Linear interpolation with proper error propagation
- **Averaging**: Inverse-variance weighted averaging in overlap regions

### Data Processing Pipeline
1. File format detection and header parsing
2. Q-range validation and optional dataset reordering
3. Overlap region identification
4. Scaling factor calculation (if enabled)
5. Merge point optimization
6. Data interpolation and averaging
7. Results validation and export

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{waldow2025saxs,
  author = {Waldow, Dean},
  title = {SAXS Data Merging Tool},
  version = {1.0.3},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.16734022},
  url = {https://doi.org/10.5281/zenodo.16734022}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Dean Waldow**
- Email: waldowda@plu.edu
- Institution: Pacific Lutheran University
- GitHub: [@waldowda](https://github.com/waldowda)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Development Setup
```bash
git clone https://github.com/waldowda/saxs-data-merging-tool.git
cd saxs-data-merging-tool
pip install -r requirements.txt
```

## Changelog

### Version 1.0.3 (Current)
- Added Q-range validation with automatic dataset swapping
- Improved data loading workflow with user feedback
- Enhanced error handling for edge cases

### Version 1.0.2
- Added error data support and uncertainty propagation
- 3-column data support (Q, I, σ)
- Error-weighted scaling and merging
- Optional error bar plotting

### Version 1.0.1
- Fixed syntax errors and plotting issues
- Full Q-range plotting implementation
- Xenocs header detection and GUI logging

### Version 1.0.0
- Initial release with GUI interface
- Automatic merge region detection
- Optional scaling functionality
- Interactive plotting and preview

## Support

For questions, bug reports, or feature requests:
- Open an issue on GitHub
- Email: waldowda@plu.edu
- Check the [documentation](https://github.com/waldowda/saxs-data-merging-tool/wiki)

---

**Keywords**: SAXS, Small Angle X-ray Scattering, Data Merging, Python, PyQt5, Scientific Computing, X-ray Analysis