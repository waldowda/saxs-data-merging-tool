"""
SAXS Data Merging Tool - PyQt5 GUI Version
===========================================

Version: 1.0.4
Created: 2025
Last Modified: August 2025

Description:
    A graphical tool for merging Small Angle X-ray Scattering (SAXS) datasets 
    with overlapping Q ranges. Features automatic merge region detection, 
    optional intensity scaling, background subtraction, and interactive visualization.

Author: Dean Waldow
Email: waldowda@plu.edu
Institution: Pacific Lutheran University

License: MIT License
    Copyright (c) 2025 Dean Waldow
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

Dependencies:
    - PyQt5 >= 5.12
    - matplotlib >= 3.3
    - numpy >= 1.19
    - pandas >= 1.2
    - scipy >= 1.6

Usage:
    python saxs_merge_gui.py

Features:
    - Interactive file selection for sample and background data
    - Background subtraction with optional scaling
    - Automatic merge region optimization
    - Optional intensity scaling
    - Real-time preview of merge regions
    - Comprehensive data export with metadata
    - Multi-threaded processing

Citation:
    If you use this tool in your research, please cite:
    Waldow, D. (2025). SAXS Data Merging Tool (Version 1.0.4) [Computer software]. 
    Zenodo. https://doi.org/10.5281/zenodo.16734022

Changelog:
    v1.0.4 - Added background subtraction capabilities
           - Optional background scaling with multiple methods
           - Enhanced GUI with background file selection
           - Improved error propagation through background subtraction
           - Added keyboard shortcuts (Ctrl+O, Ctrl+S, Ctrl+N)
    v1.0.3 - Added Q-range validation with automatic dataset swapping
           - Detects when datasets are loaded in wrong order (Dataset 1 > Dataset 2)
           - User-friendly dialog with automatic swap option
           - Improved data loading workflow and user feedback
    v1.0.2 - Added error data support and uncertainty propagation
           - 3-column data support (Q, I, σ)
           - Error-weighted scaling and merging
           - Optional error bar plotting
           - Proper uncertainty propagation in overlap regions
    v1.0.1 - Fixed syntax errors, full Q-range plotting, Xenocs header GUI logging
    v1.0.0 - Initial release with GUI interface
           - Automatic merge region detection
           - Optional scaling functionality
           - Interactive plotting and preview
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import os
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                           QWidget, QPushButton, QLabel, QLineEdit, QCheckBox, 
                           QFileDialog, QTextEdit, QGroupBox, QGridLayout, 
                           QMessageBox, QProgressBar, QDoubleSpinBox, QSpacerItem,
                           QSizePolicy, QRadioButton, QButtonGroup)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QPixmap, QKeySequence
from PyQt5.QtWidgets import QShortcut
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class BackgroundProcessor:
    """Handles background subtraction and scaling operations"""
    
    @staticmethod
    def load_background_data(filename, delimiter=None):
        """Load background data - same format as sample data"""
        return SAXSProcessor.load_saxs_data(filename, delimiter)
    
    @staticmethod
    def interpolate_background_to_sample(q_sample, q_bg, i_bg, e_bg=None):
        """Interpolate background data to sample Q points"""
        # Create interpolation function for intensity
        interp_i = interp1d(q_bg, i_bg, kind='linear', bounds_error=False, 
                          fill_value='extrapolate')
        i_bg_interp = interp_i(q_sample)
        
        # Interpolate errors if available
        if e_bg is not None:
            interp_e = interp1d(q_bg, e_bg, kind='linear', bounds_error=False, 
                              fill_value='extrapolate')
            e_bg_interp = interp_e(q_sample)
        else:
            e_bg_interp = np.sqrt(np.abs(i_bg_interp))
        
        return i_bg_interp, e_bg_interp
    
    @staticmethod
    def calculate_scaling_factor(q_sample, i_sample, q_bg, i_bg, method='auto_optimize'):
        """Calculate optimal background scaling factor"""
        if method == 'none':
            return 1.0
        
        elif method == 'auto_optimize':
            # Use high-Q region (top 20%) for optimization
            q_max = min(np.max(q_sample), np.max(q_bg))
            q_min_for_scaling = q_max * 0.8
            
            # Get data in scaling region
            mask_sample = (q_sample >= q_min_for_scaling) & (q_sample <= q_max)
            
            if np.sum(mask_sample) < 5:  # Need minimum points for optimization
                return 1.0
            
            q_scaling = q_sample[mask_sample]
            i_sample_scaling = i_sample[mask_sample]
            
            # Interpolate background to sample Q points in scaling region
            i_bg_interp, _ = BackgroundProcessor.interpolate_background_to_sample(
                q_scaling, q_bg, i_bg)
            
            # Find scaling factor that minimizes residuals
            def objective(scale_factor):
                scaled_bg = i_bg_interp * scale_factor
                residuals = np.abs(i_sample_scaling - scaled_bg)
                return np.mean(residuals)
            
            try:
                result = minimize_scalar(objective, bounds=(0.1, 10.0), method='bounded')
                return result.x
            except:
                return 1.0
        
        return 1.0
    
    @staticmethod
    def subtract_background(q_sample, i_sample, e_sample, q_bg, i_bg, e_bg=None, 
                          scaling_factor=1.0):
        """
        Subtract background from sample data with proper error propagation
        
        Returns: i_corrected, e_corrected
        """
        # Interpolate background to sample Q points
        i_bg_interp, e_bg_interp = BackgroundProcessor.interpolate_background_to_sample(
            q_sample, q_bg, i_bg, e_bg)
        
        # Apply scaling
        i_bg_scaled = i_bg_interp * scaling_factor
        e_bg_scaled = e_bg_interp * scaling_factor
        
        # Subtract background
        i_corrected = i_sample - i_bg_scaled
        
        # Propagate errors: σ² = σ_sample² + (scale × σ_bg)²
        if e_sample is not None:
            e_corrected = np.sqrt(e_sample**2 + e_bg_scaled**2)
        else:
            e_corrected = np.sqrt(np.abs(i_sample) + e_bg_scaled**2)
        
        return i_corrected, e_corrected
    
    @staticmethod
    def validate_background_subtraction(q_sample, i_sample, i_corrected, scale_factor):
        """Validate background subtraction results and return warnings"""
        warnings = []
        
        # Check scaling factor
        if scale_factor < 0.01 or scale_factor > 100:
            warnings.append(f"Unusual scaling factor: {scale_factor:.3f}")
        
        # Check for excessive negative values
        negative_fraction = np.sum(i_corrected < 0) / len(i_corrected)
        if negative_fraction > 0.1:  # More than 10% negative
            warnings.append(f"High fraction of negative intensities after subtraction: {negative_fraction*100:.1f}%")
        
        # Check for data range reduction
        intensity_reduction = (np.max(i_sample) - np.max(i_corrected)) / np.max(i_sample)
        if intensity_reduction > 0.9:  # Background is too large
            warnings.append("Background appears too large - intensity range severely reduced")
        
        return warnings

class SAXSProcessor:
    """Core SAXS processing functions - same as before but as a class"""
    
    @staticmethod
    def load_saxs_data(filename, delimiter=None):
        """
        Load SAXS data from file. Handles 2 or 3 columns (Q, I, [σ]).
        Automatically detects and skips Xenocs headers and other common formats.
        """
        try:
            # First, try to detect if this is a Xenocs file with header
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            data_start_line = 0
            xenocs_header_detected = False
            
            # Check for Xenocs header pattern
            if lines[0].strip().startswith('###########'):
                xenocs_header_detected = True
                # Find the end of the header (second line of ###########)
                for i, line in enumerate(lines[1:], 1):
                    if line.strip().startswith('###########'):
                        # Found end of header, data starts after column names
                        data_start_line = i + 2  # Skip the ########## line and column header
                        break
                # Return status message for GUI logging
                status_msg = f"Detected Xenocs header format, skipping {data_start_line} lines"
                print(status_msg)  # Also print to terminal
            
            # If Xenocs header detected, load data starting from the appropriate line
            if xenocs_header_detected:
                data = np.loadtxt(filename, delimiter=delimiter, skiprows=data_start_line)
            else:
                # Try standard loading with automatic comment detection
                data = np.loadtxt(filename, delimiter=delimiter, comments=['#', '%', ';'])
            
            # Extract Q, intensity, and optional error columns
            if data.shape[1] >= 2:
                q = data[:, 0]
                intensity = data[:, 1]
                
                # Check for error column (3rd column)
                if data.shape[1] >= 3:
                    errors = data[:, 2]
                    # Handle zero or negative errors
                    errors = np.where(errors <= 0, np.sqrt(np.abs(intensity)), errors)
                else:
                    # No error column - estimate from Poisson statistics
                    errors = np.sqrt(np.abs(intensity))
                
                return q, intensity, errors
            else:
                raise ValueError("Data file must have at least 2 columns (Q, Intensity)")
                
        except Exception as e:
            # Fallback to pandas approach
            try:
                print(f"Standard loading failed ({str(e)}), trying pandas approach...")
                
                # Try to find where data starts by looking for numeric lines
                with open(filename, 'r') as f:
                    lines = f.readlines()
                
                data_start = 0
                for i, line in enumerate(lines):
                    line = line.strip()
                    if line and not line.startswith('#') and not line.startswith('%'):
                        # Try to parse as numbers
                        try:
                            parts = line.split()
                            float(parts[0])  # Try to convert first part to float
                            float(parts[1])  # Try to convert second part to float
                            data_start = i
                            break
                        except (ValueError, IndexError):
                            continue
                
                df = pd.read_csv(filename, delimiter=delimiter, header=None, skiprows=data_start)
                
                q = df.iloc[:, 0].values
                intensity = df.iloc[:, 1].values
                
                # Check for error column
                if df.shape[1] >= 3:
                    errors = df.iloc[:, 2].values
                    errors = np.where(errors <= 0, np.sqrt(np.abs(intensity)), errors)
                else:
                    errors = np.sqrt(np.abs(intensity))
                
                return q, intensity, errors
                
            except Exception as e2:
                raise ValueError(f"Could not load data file: {str(e2)}")

    @staticmethod
    def find_overlap_range(q1, q2):
        q_min_overlap = max(np.min(q1), np.min(q2))
        q_max_overlap = min(np.max(q1), np.max(q2))
        
        if q_min_overlap >= q_max_overlap:
            raise ValueError("No overlap between datasets!")
        
        return q_min_overlap, q_max_overlap

    @staticmethod
    def calculate_scaling_factor(q1, i1, q2, i2, q_overlap_start, q_overlap_end, e1=None, e2=None):
        mask1 = (q1 >= q_overlap_start) & (q1 <= q_overlap_end)
        mask2 = (q2 >= q_overlap_start) & (q2 <= q_overlap_end)
        
        q1_overlap = q1[mask1]
        i1_overlap = i1[mask1]
        q2_overlap = q2[mask2]
        i2_overlap = i2[mask2]
        
        if len(q2_overlap) < 2:
            return 1.0
        
        interp_func = interp1d(q2_overlap, i2_overlap, kind='linear', 
                              bounds_error=False, fill_value='extrapolate')
        i2_interp = interp_func(q1_overlap)
        
        # Use error-weighted scaling if errors are provided
        if e1 is not None and e2 is not None:
            e1_overlap = e1[mask1]
            e2_overlap = e2[mask2]
            
            # Interpolate errors for dataset 2
            interp_err_func = interp1d(q2_overlap, e2_overlap, kind='linear', 
                                     bounds_error=False, fill_value='extrapolate')
            e2_interp = interp_err_func(q1_overlap)
            
            # Weight by inverse variance for more accurate scaling
            weights = 1.0 / (e1_overlap**2 + e2_interp**2)
            valid_mask = ~np.isnan(i2_interp) & (i2_interp > 0) & (i1_overlap > 0) & ~np.isnan(weights)
            
            if np.sum(valid_mask) == 0:
                return 1.0
            
            # Weighted least squares scaling factor
            weights_valid = weights[valid_mask]
            i1_valid = i1_overlap[valid_mask]
            i2_valid = i2_interp[valid_mask]
            
            scaling_factor = np.sum(weights_valid * i1_valid * i2_valid) / np.sum(weights_valid * i2_valid**2)
        else:
            # Original method without error weighting
            valid_mask = ~np.isnan(i2_interp) & (i2_interp > 0) & (i1_overlap > 0)
            if np.sum(valid_mask) == 0:
                return 1.0
            
            scaling_factor = np.median(i1_overlap[valid_mask] / i2_interp[valid_mask])
        
        return scaling_factor

    @staticmethod
    def find_best_merge_point(q1, i1, q2, i2, q_overlap_start, q_overlap_end, num_points=50):
        def merge_quality(q_merge):
            scale_factor = SAXSProcessor.calculate_scaling_factor(q1, i1, q2, i2, q_overlap_start, q_merge)
            
            window = (q_overlap_end - q_overlap_start) * 0.1
            
            mask1 = (q1 >= q_merge - window) & (q1 <= q_merge + window)
            mask2 = (q2 >= q_merge - window) & (q2 <= q_merge + window)
            
            if np.sum(mask1) < 2 or np.sum(mask2) < 2:
                return 1e6
            
            q1_local = q1[mask1]
            i1_local = i1[mask1]
            q2_local = q2[mask2]
            i2_local = i2[mask2] * scale_factor
            
            q_common = np.linspace(max(np.min(q1_local), np.min(q2_local)),
                                  min(np.max(q1_local), np.max(q2_local)), 20)
            
            if len(q_common) < 2:
                return 1e6
            
            try:
                interp1 = interp1d(q1_local, i1_local, kind='linear', bounds_error=False)
                interp2 = interp1d(q2_local, i2_local, kind='linear', bounds_error=False)
                
                i1_common = interp1(q_common)
                i2_common = interp2(q_common)
                
                valid_mask = ~(np.isnan(i1_common) | np.isnan(i2_common))
                if np.sum(valid_mask) < 2:
                    return 1e6
                
                diff = np.sqrt(np.mean((i1_common[valid_mask] - i2_common[valid_mask])**2))
                return diff
            except:
                return 1e6
        
        result = minimize_scalar(merge_quality, bounds=(q_overlap_start, q_overlap_end), method='bounded')
        return result.x

    @staticmethod
    def merge_saxs_data(q1, i1, q2, i2, q_merge_start, q_merge_end, apply_scaling=True, e1=None, e2=None):
        scale_factor = 1.0
        i2_scaled = i2.copy()
        e2_scaled = e2.copy() if e2 is not None else None
        
        if apply_scaling:
            scale_factor = SAXSProcessor.calculate_scaling_factor(q1, i1, q2, i2, q_merge_start, q_merge_end, e1, e2)
            i2_scaled = i2 * scale_factor
            if e2 is not None:
                e2_scaled = e2 * scale_factor  # Scale errors proportionally
        
        q_merged = []
        i_merged = []
        e_merged = []
        
        # Add dataset 1 up to merge start
        mask1_before = q1 < q_merge_start
        q_merged.extend(q1[mask1_before])
        i_merged.extend(i1[mask1_before])
        if e1 is not None:
            e_merged.extend(e1[mask1_before])
        else:
            e_merged.extend(np.sqrt(np.abs(i1[mask1_before])))
        
        # Create overlap region with proper error propagation
        q_overlap = np.linspace(q_merge_start, q_merge_end, 100)
        
        interp1 = interp1d(q1, i1, kind='linear', bounds_error=False, fill_value=np.nan)
        interp2 = interp1d(q2, i2_scaled, kind='linear', bounds_error=False, fill_value=np.nan)
        
        # Interpolate errors
        if e1 is not None:
            interp_e1 = interp1d(q1, e1, kind='linear', bounds_error=False, fill_value=np.nan)
        else:
            interp_e1 = lambda q: np.sqrt(np.abs(interp1(q)))
            
        if e2_scaled is not None:
            interp_e2 = interp1d(q2, e2_scaled, kind='linear', bounds_error=False, fill_value=np.nan)
        else:
            interp_e2 = lambda q: np.sqrt(np.abs(interp2(q)))
        
        i1_overlap = interp1(q_overlap)
        i2_overlap = interp2(q_overlap)
        e1_overlap = interp_e1(q_overlap)
        e2_overlap = interp_e2(q_overlap)
        
        # Weighted averaging in overlap region
        for i, q_val in enumerate(q_overlap):
            if not np.isnan(i1_overlap[i]) and not np.isnan(i2_overlap[i]):
                # Error-weighted average
                if not np.isnan(e1_overlap[i]) and not np.isnan(e2_overlap[i]) and e1_overlap[i] > 0 and e2_overlap[i] > 0:
                    w1 = 1.0 / e1_overlap[i]**2
                    w2 = 1.0 / e2_overlap[i]**2
                    w_total = w1 + w2
                    
                    i_avg = (w1 * i1_overlap[i] + w2 * i2_overlap[i]) / w_total
                    e_avg = np.sqrt(1.0 / w_total)  # Propagated uncertainty
                else:
                    # Simple average if errors unavailable
                    i_avg = (i1_overlap[i] + i2_overlap[i]) / 2
                    e_avg = np.sqrt(e1_overlap[i]**2 + e2_overlap[i]**2) / 2
                
                q_merged.append(q_val)
                i_merged.append(i_avg)
                e_merged.append(e_avg)
                
            elif not np.isnan(i1_overlap[i]):
                q_merged.append(q_val)
                i_merged.append(i1_overlap[i])
                e_merged.append(e1_overlap[i] if not np.isnan(e1_overlap[i]) else np.sqrt(abs(i1_overlap[i])))
            elif not np.isnan(i2_overlap[i]):
                q_merged.append(q_val)
                i_merged.append(i2_overlap[i])
                e_merged.append(e2_overlap[i] if not np.isnan(e2_overlap[i]) else np.sqrt(abs(i2_overlap[i])))
        
        # Add dataset 2 after merge end
        mask2_after = q2 > q_merge_end
        q_merged.extend(q2[mask2_after])
        i_merged.extend(i2_scaled[mask2_after])
        if e2_scaled is not None:
            e_merged.extend(e2_scaled[mask2_after])
        else:
            e_merged.extend(np.sqrt(np.abs(i2_scaled[mask2_after])))
        
        # Convert to arrays and sort by Q
        q_merged = np.array(q_merged)
        i_merged = np.array(i_merged)
        e_merged = np.array(e_merged)
        
        sort_indices = np.argsort(q_merged)
        q_merged = q_merged[sort_indices]
        i_merged = i_merged[sort_indices]
        e_merged = e_merged[sort_indices]
        
        return q_merged, i_merged, e_merged, scale_factor

class MergeWorker(QThread):
    """Worker thread for processing to avoid GUI freezing"""
    finished = pyqtSignal(object, object, object, float)  # Added error data
    progress = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, q1, i1, q2, i2, q_merge_start, q_merge_end, apply_scaling, e1=None, e2=None):
        super().__init__()
        self.q1, self.i1 = q1, i1
        self.q2, self.i2 = q2, i2
        self.e1, self.e2 = e1, e2
        self.q_merge_start = q_merge_start
        self.q_merge_end = q_merge_end
        self.apply_scaling = apply_scaling
    
    def run(self):
        try:
            self.progress.emit("Merging datasets...")
            q_merged, i_merged, e_merged, scale_factor = SAXSProcessor.merge_saxs_data(
                self.q1, self.i1, self.q2, self.i2, 
                self.q_merge_start, self.q_merge_end, self.apply_scaling, self.e1, self.e2)
            self.finished.emit(q_merged, i_merged, e_merged, scale_factor)
        except Exception as e:
            self.error.emit(str(e))

class PlotCanvas(FigureCanvas):
    """Matplotlib canvas for displaying plots in Qt"""
    
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(12, 8))
        super().__init__(self.fig)
        self.setParent(parent)
        
    def plot_background_comparison(self, q1, i1, q_bg1, i_bg1, bg1_scale,
                                 q2, i2, q_bg2, i_bg2, bg2_scale,
                                 show_errors=True, e1=None, e_bg1=None, 
                                 e2=None, e_bg2=None):
        """Plot sample and background data for comparison"""
        self.fig.clear()
        
        # Dataset 1
        ax1 = self.fig.add_subplot(2, 1, 1)
        
        if show_errors and e1 is not None and e_bg1 is not None:
            ax1.errorbar(q1, i1, yerr=e1, fmt='b-', label='Sample 1', alpha=0.8,
                        errorevery=max(1, len(q1)//30), capsize=2)
            ax1.errorbar(q_bg1, i_bg1, yerr=e_bg1, fmt='r--', label='Background 1', alpha=0.6,
                        errorevery=max(1, len(q_bg1)//30), capsize=2)
            if bg1_scale != 1.0:
                ax1.errorbar(q_bg1, i_bg1 * bg1_scale, yerr=e_bg1 * bg1_scale, 
                           fmt='r:', label=f'Background 1 (scaled {bg1_scale:.3f})', alpha=0.8,
                           errorevery=max(1, len(q_bg1)//30), capsize=2)
        else:
            ax1.loglog(q1, i1, 'b-', label='Sample 1', alpha=0.8)
            ax1.loglog(q_bg1, i_bg1, 'r--', label='Background 1', alpha=0.6)
            if bg1_scale != 1.0:
                ax1.loglog(q_bg1, i_bg1 * bg1_scale, 'r:', 
                         label=f'Background 1 (scaled {bg1_scale:.3f})', alpha=0.8)
        
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Q (Å⁻¹)')
        ax1.set_ylabel('Intensity')
        ax1.set_title('Dataset 1: Sample vs Background')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Dataset 2
        ax2 = self.fig.add_subplot(2, 1, 2)
        
        if show_errors and e2 is not None and e_bg2 is not None:
            ax2.errorbar(q2, i2, yerr=e2, fmt='b-', label='Sample 2', alpha=0.8,
                        errorevery=max(1, len(q2)//30), capsize=2)
            ax2.errorbar(q_bg2, i_bg2, yerr=e_bg2, fmt='r--', label='Background 2', alpha=0.6,
                        errorevery=max(1, len(q_bg2)//30), capsize=2)
            if bg2_scale != 1.0:
                ax2.errorbar(q_bg2, i_bg2 * bg2_scale, yerr=e_bg2 * bg2_scale, 
                           fmt='r:', label=f'Background 2 (scaled {bg2_scale:.3f})', alpha=0.8,
                           errorevery=max(1, len(q_bg2)//30), capsize=2)
        else:
            ax2.loglog(q2, i2, 'b-', label='Sample 2', alpha=0.8)
            ax2.loglog(q_bg2, i_bg2, 'r--', label='Background 2', alpha=0.6)
            if bg2_scale != 1.0:
                ax2.loglog(q_bg2, i_bg2 * bg2_scale, 'r:', 
                         label=f'Background 2 (scaled {bg2_scale:.3f})', alpha=0.8)
        
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel('Q (Å⁻¹)')
        ax2.set_ylabel('Intensity')
        ax2.set_title('Dataset 2: Sample vs Background')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        
    def plot_overlap_region(self, q1, i1, q2, i2, q_overlap_start, q_overlap_end, 
                           q_merge_start=None, q_merge_end=None, scale_factor=1.0, 
                           show_errors=True, e1=None, e2=None):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        
        # Plot full datasets (FIXED: Show complete Q range)
        if show_errors and e1 is not None and e2 is not None:
            # Plot with error bars
            ax.errorbar(q1, i1, yerr=e1, fmt='b-', label='Dataset 1', linewidth=2, alpha=0.7, 
                       errorevery=max(1, len(q1)//50), capsize=3)
            ax.errorbar(q2, i2, yerr=e2, fmt='r-', label='Dataset 2 (original)', alpha=0.5,
                       errorevery=max(1, len(q2)//50), capsize=3)
            ax.errorbar(q2, i2 * scale_factor, yerr=e2 * scale_factor, fmt='r--', 
                       label='Dataset 2 (scaled)', linewidth=2, 
                       errorevery=max(1, len(q2)//50), capsize=3)
        else:
            # Plot without error bars
            ax.loglog(q1, i1, 'b-', label='Dataset 1', linewidth=2, alpha=0.7)
            ax.loglog(q2, i2, 'r-', label='Dataset 2 (original)', alpha=0.5)
            ax.loglog(q2, i2 * scale_factor, 'r--', label='Dataset 2 (scaled)', linewidth=2)
        
        # Set log scale
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        # Highlight overlap region
        ax.axvspan(q_overlap_start, q_overlap_end, alpha=0.2, color='gray', label='Available overlap')
        
        # Highlight proposed merge region if provided
        if q_merge_start is not None and q_merge_end is not None:
            ax.axvspan(q_merge_start, q_merge_end, alpha=0.3, color='green', label='Proposed merge region')
            ax.axvline(q_merge_start, color='g', linestyle=':', linewidth=2)
            ax.axvline(q_merge_end, color='g', linestyle=':', linewidth=2)
        
        # Set axis limits to show full data range (FIXED: Ensure full range visible)
        q_min = min(np.min(q1), np.min(q2))
        q_max = max(np.max(q1), np.max(q2))
        ax.set_xlim(q_min * 0.8, q_max * 1.2)  # Add some padding
        
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('Intensity (log scale)')
        ax.set_title('Full Dataset View - Assess Merge Quality')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        self.draw()
    
    def plot_final_results(self, q1, i1, q2, i2, q_merged, i_merged, 
                          q_merge_start, q_merge_end, scale_factor, apply_scaling,
                          show_errors=True, e1=None, e2=None, e_merged=None):
        self.fig.clear()
        
        # Original datasets
        ax1 = self.fig.add_subplot(2, 1, 1)
        
        if show_errors and e1 is not None and e2 is not None:
            # Plot with error bars
            ax1.errorbar(q1, i1, yerr=e1, fmt='b-', label='Dataset 1', alpha=0.7,
                        errorevery=max(1, len(q1)//30), capsize=2)
            ax1.errorbar(q2, i2, yerr=e2, fmt='r-', label='Dataset 2', alpha=0.7,
                        errorevery=max(1, len(q2)//30), capsize=2)
            if apply_scaling:
                ax1.errorbar(q2, i2 * scale_factor, yerr=e2 * scale_factor, fmt='r--', 
                           label='Dataset 2 (scaled)', alpha=0.7,
                           errorevery=max(1, len(q2)//30), capsize=2)
        else:
            # Plot without error bars
            ax1.loglog(q1, i1, 'b-', label='Dataset 1', alpha=0.7)
            ax1.loglog(q2, i2, 'r-', label='Dataset 2', alpha=0.7)
            if apply_scaling:
                ax1.loglog(q2, i2 * scale_factor, 'r--', label='Dataset 2 (scaled)', alpha=0.7)
        
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.axvline(q_merge_start, color='g', linestyle=':', label='Merge start')
        ax1.axvline(q_merge_end, color='g', linestyle=':', label='Merge end')
        ax1.set_xlabel('Q (Å⁻¹)')
        ax1.set_ylabel('Intensity (log scale)')
        title_text = 'Original SAXS Datasets (Log-Log Plot)'
        if apply_scaling:
            title_text += f' - Scaling factor: {scale_factor:.4f}'
        else:
            title_text += ' - No scaling applied'
        ax1.set_title(title_text)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Merged dataset
        ax2 = self.fig.add_subplot(2, 1, 2)
        
        if show_errors and e_merged is not None:
            # Plot merged data with error bars
            ax2.errorbar(q_merged, i_merged, yerr=e_merged, fmt='k-', linewidth=2, 
                        label='Merged dataset', errorevery=max(1, len(q_merged)//50), capsize=3)
            # Background datasets without error bars for clarity
            ax2.loglog(q1, i1, 'b-', alpha=0.3, label='Dataset 1')
            if apply_scaling:
                ax2.loglog(q2, i2 * scale_factor, 'r-', alpha=0.3, label='Dataset 2 (scaled)')
            else:
                ax2.loglog(q2, i2, 'r-', alpha=0.3, label='Dataset 2')
        else:
            # Plot without error bars
            ax2.loglog(q_merged, i_merged, 'k-', linewidth=2, label='Merged dataset')
            ax2.loglog(q1, i1, 'b-', alpha=0.3, label='Dataset 1')
            if apply_scaling:
                ax2.loglog(q2, i2 * scale_factor, 'r-', alpha=0.3, label='Dataset 2 (scaled)')
            else:
                ax2.loglog(q2, i2, 'r-', alpha=0.3, label='Dataset 2')
        
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.axvline(q_merge_start, color='g', linestyle=':', label='Merge start')
        ax2.axvline(q_merge_end, color='g', linestyle=':', label='Merge end')
        ax2.set_xlabel('Q (Å⁻¹)')
        ax2.set_ylabel('Intensity (log scale)')
        ax2.set_title('Merged SAXS Dataset (Log-Log Plot)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()

class SAXSMergeGUI(QMainWindow):
    """Main GUI application"""
    
    def __init__(self):
        super().__init__()
        # Sample data
        self.q1 = self.i1 = self.e1 = self.q2 = self.i2 = self.e2 = None
        # Background data
        self.q_bg1 = self.i_bg1 = self.e_bg1 = None
        self.q_bg2 = self.i_bg2 = self.e_bg2 = None
        # Overlap and merge data
        self.q_overlap_start = self.q_overlap_end = None
        self.q_merged = self.i_merged = self.e_merged = None
        self.scale_factor = 1.0
        # File paths
        self.file1 = self.file2 = ""
        self.bg_file1 = self.bg_file2 = ""
        # Background scaling factors
        self.bg1_scale = self.bg2_scale = 1.0
        
        self.initUI()
        self.setup_shortcuts()
        
    def initUI(self):
        self.setWindowTitle('SAXS Data Merging Tool v1.0.4')
        self.setGeometry(100, 100, 1400, 1000)
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # Title and info section
        title_widget = QWidget()
        title_layout = QVBoxLayout(title_widget)
        
        # Main title
        title = QLabel('SAXS Data Merging Tool')
        title.setFont(QFont('Arial', 16, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        title_layout.addWidget(title)
        
        # Version info
        version_label = QLabel('Version 1.0.4 - Now with Background Subtraction!')
        version_label.setFont(QFont('Arial', 10))
        version_label.setAlignment(Qt.AlignCenter)
        version_label.setStyleSheet("color: #0066cc; font-weight: bold;")
        title_layout.addWidget(version_label)
        
        # Contact/License info space
        info_label = QLabel(
            'Author: Dean Waldow | Email: waldowda@plu.edu\n'
            'Institution: Pacific Lutheran University | License: MIT\n'
            'DOI: 10.5281/zenodo.16734022 | Keyboard shortcuts: Ctrl+O (Open), Ctrl+S (Save), Ctrl+N (New)'
        )
        info_label.setFont(QFont('Arial', 9))
        info_label.setAlignment(Qt.AlignCenter)
        info_label.setStyleSheet("color: #666666; padding: 5px;")
        info_label.setWordWrap(True)
        title_layout.addWidget(info_label)
        
        title_layout.addSpacing(10)
        layout.addWidget(title_widget)
        
        # File selection group
        file_group = QGroupBox("File Selection")
        file_layout = QGridLayout(file_group)
        
        self.file1_label = QLabel("No file selected")
        self.file2_label = QLabel("No file selected")
        self.bg_file1_label = QLabel("No background file")
        self.bg_file2_label = QLabel("No background file")
        
        # Sample files
        file_layout.addWidget(QLabel("Dataset 1 (lower Q data):"), 0, 0)
        file_layout.addWidget(self.file1_label, 0, 1)
        file_layout.addWidget(QPushButton("Browse...", clicked=self.select_file1), 0, 2)
        
        file_layout.addWidget(QLabel("Dataset 2 (higher Q data):"), 1, 0)
        file_layout.addWidget(self.file2_label, 1, 1)
        file_layout.addWidget(QPushButton("Browse...", clicked=self.select_file2), 1, 2)
        
        # Background files
        file_layout.addWidget(QLabel("Background 1:"), 2, 0)
        file_layout.addWidget(self.bg_file1_label, 2, 1)
        bg1_layout = QHBoxLayout()
        bg1_btn = QPushButton("Browse...", clicked=self.select_bg_file1)
        clear_bg1_btn = QPushButton("Clear", clicked=self.clear_bg_file1)
        clear_bg1_btn.setMaximumWidth(60)
        bg1_layout.addWidget(bg1_btn)
        bg1_layout.addWidget(clear_bg1_btn)
        file_layout.addLayout(bg1_layout, 2, 2)
        
        file_layout.addWidget(QLabel("Background 2:"), 3, 0)
        file_layout.addWidget(self.bg_file2_label, 3, 1)
        bg2_layout = QHBoxLayout()
        bg2_btn = QPushButton("Browse...", clicked=self.select_bg_file2)
        clear_bg2_btn = QPushButton("Clear", clicked=self.clear_bg_file2)
        clear_bg2_btn.setMaximumWidth(60)
        bg2_layout.addWidget(bg2_btn)
        bg2_layout.addWidget(clear_bg2_btn)
        file_layout.addLayout(bg2_layout, 3, 2)
        
        layout.addWidget(file_group)
        
        # Background options group
        bg_options_group = QGroupBox("Background Subtraction Options")
        bg_options_layout = QVBoxLayout(bg_options_group)
        
        self.enable_bg_subtraction = QCheckBox("Enable background subtraction")
        self.enable_bg_subtraction.setChecked(False)
        self.enable_bg_subtraction.stateChanged.connect(self.on_bg_subtraction_toggle)
        bg_options_layout.addWidget(self.enable_bg_subtraction)
        
        self.same_bg_checkbox = QCheckBox("Use same background for both datasets (rare)")
        self.same_bg_checkbox.setEnabled(False)
        bg_options_layout.addWidget(self.same_bg_checkbox)
        
        self.enable_bg_scaling = QCheckBox("Apply scaling to backgrounds before subtraction")
        self.enable_bg_scaling.setChecked(False)  # Default OFF
        self.enable_bg_scaling.setEnabled(False)
        self.enable_bg_scaling.stateChanged.connect(self.on_bg_scaling_toggle)
        bg_options_layout.addWidget(self.enable_bg_scaling)
        
        # Background scaling options (initially hidden)
        self.bg_scaling_widget = QWidget()
        bg_scaling_layout = QHBoxLayout(self.bg_scaling_widget)
        bg_scaling_layout.addWidget(QLabel("Scaling method:"))
        
        self.bg_scaling_method = QButtonGroup()
        
        no_scaling = QRadioButton("No scaling")
        no_scaling.setChecked(True)
        self.bg_scaling_method.addButton(no_scaling, 0)
        bg_scaling_layout.addWidget(no_scaling)
        
        auto_scaling = QRadioButton("Auto optimize")
        self.bg_scaling_method.addButton(auto_scaling, 1)
        bg_scaling_layout.addWidget(auto_scaling)
        
        bg_scaling_layout.addWidget(QLabel("Manual scales:"))
        bg_scaling_layout.addWidget(QLabel("BG1:"))
        self.bg1_scale_spin = QDoubleSpinBox()
        self.bg1_scale_spin.setRange(0.001, 1000.0)
        self.bg1_scale_spin.setValue(1.0)
        self.bg1_scale_spin.setDecimals(4)
        self.bg1_scale_spin.setMaximumWidth(80)
        bg_scaling_layout.addWidget(self.bg1_scale_spin)
        
        bg_scaling_layout.addWidget(QLabel("BG2:"))
        self.bg2_scale_spin = QDoubleSpinBox()
        self.bg2_scale_spin.setRange(0.001, 1000.0)
        self.bg2_scale_spin.setValue(1.0)
        self.bg2_scale_spin.setDecimals(4)
        self.bg2_scale_spin.setMaximumWidth(80)
        bg_scaling_layout.addWidget(self.bg2_scale_spin)
        
        bg_scaling_layout.addStretch()
        
        self.bg_scaling_widget.setVisible(False)
        bg_options_layout.addWidget(self.bg_scaling_widget)
        
        layout.addWidget(bg_options_group)
        
        # General options group
        options_group = QGroupBox("General Options")
        options_layout = QVBoxLayout(options_group)
        
        self.scaling_checkbox = QCheckBox("Apply intensity scaling to match datasets")
        self.scaling_checkbox.setChecked(True)
        options_layout.addWidget(self.scaling_checkbox)
        
        self.errorbar_checkbox = QCheckBox("Show error bars in plots")
        self.errorbar_checkbox.setChecked(True)
        self.errorbar_checkbox.stateChanged.connect(self.on_errorbar_toggle)
        options_layout.addWidget(self.errorbar_checkbox)
        
        self.show_bg_plots = QCheckBox("Show background data in plots")
        self.show_bg_plots.setChecked(False)
        self.show_bg_plots.setEnabled(False)
        options_layout.addWidget(self.show_bg_plots)
        
        layout.addWidget(options_group)
        
        # Merge region group
        merge_group = QGroupBox("Merge Region")
        merge_layout = QGridLayout(merge_group)
        
        merge_layout.addWidget(QLabel("Start Q:"), 0, 0)
        self.q_start_spin = QDoubleSpinBox()
        self.q_start_spin.setDecimals(6)
        self.q_start_spin.setRange(0.0, 10.0)
        self.q_start_spin.setSingleStep(0.001)
        merge_layout.addWidget(self.q_start_spin, 0, 1)
        
        merge_layout.addWidget(QLabel("End Q:"), 0, 2)
        self.q_end_spin = QDoubleSpinBox()
        self.q_end_spin.setDecimals(6)
        self.q_end_spin.setRange(0.0, 10.0)
        self.q_end_spin.setSingleStep(0.001)
        merge_layout.addWidget(self.q_end_spin, 0, 3)
        
        self.auto_merge_btn = QPushButton("Find Auto Merge Region", clicked=self.find_auto_merge)
        self.preview_btn = QPushButton("Preview Merge Region", clicked=self.preview_merge)
        
        merge_layout.addWidget(self.auto_merge_btn, 1, 0, 1, 2)
        merge_layout.addWidget(self.preview_btn, 1, 2, 1, 2)
        
        layout.addWidget(merge_group)
        
        # Plot canvas
        self.plot_canvas = PlotCanvas()
        layout.addWidget(self.plot_canvas)
        
        # Action buttons
        button_layout = QHBoxLayout()
        
        self.new_analysis_btn = QPushButton("New Analysis (Ctrl+N)", clicked=self.reset_analysis)
        self.merge_btn = QPushButton("Merge Datasets", clicked=self.merge_datasets)
        self.merge_btn.setEnabled(False)
        
        self.save_btn = QPushButton("Save Results (Ctrl+S)", clicked=self.save_results)
        self.save_btn.setEnabled(False)
        
        # Background specific buttons
        self.preview_bg_btn = QPushButton("Preview Background Subtraction", clicked=self.preview_background)
        self.preview_bg_btn.setEnabled(False)
        
        button_layout.addWidget(self.new_analysis_btn)
        button_layout.addWidget(self.preview_bg_btn)
        button_layout.addWidget(self.merge_btn)
        button_layout.addWidget(self.save_btn)
        button_layout.addStretch()
        
        layout.addLayout(button_layout)
        
        # Status and progress
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(120)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
    
    def setup_shortcuts(self):
        """Setup keyboard shortcuts"""
        # Ctrl+O for open file
        open_shortcut = QShortcut(QKeySequence("Ctrl+O"), self)
        open_shortcut.activated.connect(self.select_file1)  # Default to first file
        
        # Ctrl+S for save
        save_shortcut = QShortcut(QKeySequence("Ctrl+S"), self)
        save_shortcut.activated.connect(self.save_results)
        
        # Ctrl+N for new analysis
        new_shortcut = QShortcut(QKeySequence("Ctrl+N"), self)
        new_shortcut.activated.connect(self.reset_analysis)
    
    def log_message(self, message):
        """Add message to status text"""
        self.status_text.append(message)
        QApplication.processEvents()
    
    def on_bg_subtraction_toggle(self):
        """Handle background subtraction enable/disable"""
        enabled = self.enable_bg_subtraction.isChecked()
        self.same_bg_checkbox.setEnabled(enabled)
        self.enable_bg_scaling.setEnabled(enabled)
        self.show_bg_plots.setEnabled(enabled)
        
        # Enable preview if at least one background file is loaded
        bg_files_available = bool(self.bg_file1) or bool(self.bg_file2)
        self.preview_bg_btn.setEnabled(enabled and bg_files_available)
        
        if enabled:
            self.log_message("Background subtraction enabled")
        else:
            self.log_message("Background subtraction disabled")
            
    def on_bg_scaling_toggle(self):
        """Handle background scaling enable/disable"""
        enabled = self.enable_bg_scaling.isChecked()
        self.bg_scaling_widget.setVisible(enabled)
        
        if enabled:
            self.log_message("Background scaling enabled")
        else:
            self.log_message("Background scaling disabled")
    
    def reset_analysis(self):
        """Reset all analysis data and UI state for new analysis"""
        # Stop any running worker
        if hasattr(self, 'worker') and self.worker.isRunning():
            self.worker.terminate()
            self.worker.wait()
        
        # Clear sample data
        self.q1 = self.i1 = self.e1 = self.q2 = self.i2 = self.e2 = None
        # Clear background data
        self.q_bg1 = self.i_bg1 = self.e_bg1 = None
        self.q_bg2 = self.i_bg2 = self.e_bg2 = None
        self.q_overlap_start = self.q_overlap_end = None
        self.q_merged = self.i_merged = self.e_merged = None
        self.scale_factor = 1.0
        self.bg1_scale = self.bg2_scale = 1.0
        
        # Reset file paths and labels
        self.file1 = self.file2 = ""
        self.bg_file1 = self.bg_file2 = ""
        self.file1_label.setText("No file selected")
        self.file2_label.setText("No file selected")
        self.bg_file1_label.setText("No background file")
        self.bg_file2_label.setText("No background file")
        
        # Reset UI state
        self.q_start_spin.setRange(0.0, 10.0)
        self.q_end_spin.setRange(0.0, 10.0)
        self.q_start_spin.setValue(0.0)
        self.q_end_spin.setValue(0.0)
        
        # Reset checkboxes to defaults
        self.scaling_checkbox.setChecked(True)
        self.enable_bg_subtraction.setChecked(False)
        self.enable_bg_scaling.setChecked(False)
        self.same_bg_checkbox.setChecked(False)
        self.bg1_scale_spin.setValue(1.0)
        self.bg2_scale_spin.setValue(1.0)
        
        # Disable controls
        self.auto_merge_btn.setEnabled(False)
        self.preview_btn.setEnabled(False)
        self.merge_btn.setEnabled(False)
        self.save_btn.setEnabled(False)
        self.preview_bg_btn.setEnabled(False)
        
        # Clear plot
        self.plot_canvas.fig.clear()
        self.plot_canvas.draw()
        
        # Hide progress bar
        self.progress_bar.setVisible(False)
        
        # Clear status log
        self.status_text.clear()
        self.log_message("--- New Analysis Started ---")
        self.log_message("v1.0.4: Background subtraction feature available!")
    
    def on_errorbar_toggle(self):
        """Handle error bar checkbox toggle - refresh current plot"""
        if hasattr(self, 'q_merged') and self.q_merged is not None:
            # Refresh final results plot
            self.plot_canvas.plot_final_results(
                self.q1, self.i1, self.q2, self.i2, self.q_merged, self.i_merged,
                self.q_start_spin.value(), self.q_end_spin.value(), 
                self.scale_factor, self.scaling_checkbox.isChecked(),
                self.errorbar_checkbox.isChecked(), self.e1, self.e2, self.e_merged)
        elif hasattr(self, 'q1') and self.q1 is not None and hasattr(self, 'q2') and self.q2 is not None:
            # Refresh preview plot
            self.preview_merge()
        
    def select_file1(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select first SAXS dataset", "",
            "Data files (*.dat *.txt *.csv *.xy *.chi);;All files (*.*)")
        
        if file_path:
            self.file1 = file_path
            self.file1_label.setText(os.path.basename(file_path))
            self.log_message(f"Selected dataset 1: {os.path.basename(file_path)}")
            self.load_and_analyze()
    
    def select_file2(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select second SAXS dataset", "",
            "Data files (*.dat *.txt *.csv *.xy *.chi);;All files (*.*)")
        
        if file_path:
            self.file2 = file_path
            self.file2_label.setText(os.path.basename(file_path))
            self.log_message(f"Selected dataset 2: {os.path.basename(file_path)}")
            self.load_and_analyze()
    
    def select_bg_file1(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select background for dataset 1", "",
            "Data files (*.dat *.txt *.csv *.xy *.chi);;All files (*.*)")
        
        if file_path:
            self.bg_file1 = file_path
            self.bg_file1_label.setText(os.path.basename(file_path))
            self.log_message(f"Selected background 1: {os.path.basename(file_path)}")
            self.load_background_data()
    
    def select_bg_file2(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select background for dataset 2", "",
            "Data files (*.dat *.txt *.csv *.xy *.chi);;All files (*.*)")
        
        if file_path:
            self.bg_file2 = file_path
            self.bg_file2_label.setText(os.path.basename(file_path))
            self.log_message(f"Selected background 2: {os.path.basename(file_path)}")
            self.load_background_data()
    
    def clear_bg_file1(self):
        self.bg_file1 = ""
        self.bg_file1_label.setText("No background file")
        self.q_bg1 = self.i_bg1 = self.e_bg1 = None
        self.log_message("Cleared background 1")
        self.preview_bg_btn.setEnabled(False)
    
    def clear_bg_file2(self):
        self.bg_file2 = ""
        self.bg_file2_label.setText("No background file")
        self.q_bg2 = self.i_bg2 = self.e_bg2 = None
        self.log_message("Cleared background 2")
        self.preview_bg_btn.setEnabled(False)
    
    def load_background_data(self):
        """Load background data files"""
        try:
            if self.bg_file1:
                self.q_bg1, self.i_bg1, self.e_bg1 = BackgroundProcessor.load_background_data(self.bg_file1)
                self.log_message(f"Background 1: Q range {self.q_bg1.min():.4f} to {self.q_bg1.max():.4f} ({len(self.q_bg1)} points)")
                
                # Check for error data
                if self.e_bg1 is not None and np.any(self.e_bg1 != np.sqrt(np.abs(self.i_bg1))):
                    self.log_message("Background 1: Error column detected")
                else:
                    self.log_message("Background 1: Using Poisson error estimates")
                    
            if self.bg_file2:
                self.q_bg2, self.i_bg2, self.e_bg2 = BackgroundProcessor.load_background_data(self.bg_file2)
                self.log_message(f"Background 2: Q range {self.q_bg2.min():.4f} to {self.q_bg2.max():.4f} ({len(self.q_bg2)} points)")
                
                # Check for error data
                if self.e_bg2 is not None and np.any(self.e_bg2 != np.sqrt(np.abs(self.i_bg2))):
                    self.log_message("Background 2: Error column detected")
                else:
                    self.log_message("Background 2: Using Poisson error estimates")
            
            # Enable preview button if any background files loaded
            bg_files_available = bool(self.bg_file1) or bool(self.bg_file2)
            if bg_files_available and self.enable_bg_subtraction.isChecked():
                self.preview_bg_btn.setEnabled(True)
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading background data: {str(e)}")
    
    def load_and_analyze(self):
        """Load data and analyze overlap if both files selected"""
        if not (self.file1 and self.file2):
            return
            
        try:
            # Load sample data
            self.log_message("Loading SAXS datasets...")
            self.q1, self.i1, self.e1 = SAXSProcessor.load_saxs_data(self.file1)
            
            # Check for Xenocs headers
            if self.file1:
                with open(self.file1, 'r') as f:
                    first_line = f.readline().strip()
                if first_line.startswith('###########'):
                    self.log_message("Detected Xenocs header in dataset 1 - automatically skipped")
            
            self.q2, self.i2, self.e2 = SAXSProcessor.load_saxs_data(self.file2)
            
            if self.file2:
                with open(self.file2, 'r') as f:
                    first_line = f.readline().strip()
                if first_line.startswith('###########'):
                    self.log_message("Detected Xenocs header in dataset 2 - automatically skipped")
            
            self.log_message(f"Dataset 1: Q range {self.q1.min():.4f} to {self.q1.max():.4f} ({len(self.q1)} points)")
            self.log_message(f"Dataset 2: Q range {self.q2.min():.4f} to {self.q2.max():.4f} ({len(self.q2)} points)")
            
            # Validate Q range ordering
            q1_mean = np.mean(self.q1)
            q2_mean = np.mean(self.q2)
            
            if q1_mean > q2_mean:
                reply = QMessageBox.question(
                    self, "Q Range Warning", 
                    f"Dataset 1 has higher average Q ({q1_mean:.4f}) than Dataset 2 ({q2_mean:.4f}).\n\n"
                    "For optimal merging, Dataset 1 should contain lower Q data and Dataset 2 should contain higher Q data.\n\n"
                    "Would you like to swap the datasets automatically?",
                    QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
                    QMessageBox.Yes
                )
                
                if reply == QMessageBox.Yes:
                    # Swap the datasets
                    self.file1, self.file2 = self.file2, self.file1
                    self.q1, self.i1, self.e1, self.q2, self.i2, self.e2 = self.q2, self.i2, self.e2, self.q1, self.i1, self.e1
                    
                    # Update labels
                    self.file1_label.setText(os.path.basename(self.file1))
                    self.file2_label.setText(os.path.basename(self.file2))
                    
                    self.log_message("Datasets swapped automatically.")
                    self.log_message(f"Dataset 1 (lower Q): Q range {self.q1.min():.4f} to {self.q1.max():.4f}")
                    self.log_message(f"Dataset 2 (higher Q): Q range {self.q2.min():.4f} to {self.q2.max():.4f}")
                    
                elif reply == QMessageBox.Cancel:
                    self.log_message("Data loading cancelled.")
                    return
                else:
                    self.log_message("Proceeding with current dataset order (may not be optimal).")
            
            # Check error data
            if self.e1 is not None and np.any(self.e1 != np.sqrt(np.abs(self.i1))):
                self.log_message("Dataset 1: Error column detected and loaded")
            else:
                self.log_message("Dataset 1: No error column - using Poisson estimates")
                
            if self.e2 is not None and np.any(self.e2 != np.sqrt(np.abs(self.i2))):
                self.log_message("Dataset 2: Error column detected and loaded")
            else:
                self.log_message("Dataset 2: No error column - using Poisson estimates")
            
            # Find overlap
            self.q_overlap_start, self.q_overlap_end = SAXSProcessor.find_overlap_range(self.q1, self.q2)
            self.log_message(f"Overlap range: {self.q_overlap_start:.4f} to {self.q_overlap_end:.4f}")
            
            # Enable controls
            self.auto_merge_btn.setEnabled(True)
            self.preview_btn.setEnabled(True)
            self.merge_btn.setEnabled(True)
            
            # Set initial merge range
            self.q_start_spin.setRange(self.q_overlap_start, self.q_overlap_end)
            self.q_end_spin.setRange(self.q_overlap_start, self.q_overlap_end)
            
            # Enable background preview if any backgrounds loaded
            bg_files_available = bool(self.bg_file1) or bool(self.bg_file2)
            if bg_files_available and self.enable_bg_subtraction.isChecked():
                self.preview_bg_btn.setEnabled(True)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading data: {str(e)}")
    
    def apply_background_subtraction(self):
        """Apply background subtraction to sample data"""
        if not self.enable_bg_subtraction.isChecked():
            return self.q1, self.i1, self.e1, self.q2, self.i2, self.e2
        
        # Check if any background data is available
        if not (self.q_bg1 is not None or self.q_bg2 is not None):
            self.log_message("Warning: Background subtraction enabled but no background data loaded")
            return self.q1, self.i1, self.e1, self.q2, self.i2, self.e2
        
        try:
            # Initialize corrected data with original data
            q1_corrected, i1_corrected, e1_corrected = self.q1, self.i1, self.e1
            q2_corrected, i2_corrected, e2_corrected = self.q2, self.i2, self.e2
            
            # Process dataset 1 if background available
            if self.q_bg1 is not None:
                # Calculate scaling factor if enabled
                if self.enable_bg_scaling.isChecked():
                    scaling_method_id = self.bg_scaling_method.checkedId()
                    
                    if scaling_method_id == 0:  # No scaling
                        self.bg1_scale = 1.0
                    elif scaling_method_id == 1:  # Auto optimize
                        self.bg1_scale = BackgroundProcessor.calculate_scaling_factor(
                            self.q1, self.i1, self.q_bg1, self.i_bg1, 'auto_optimize')
                    else:  # Manual scaling (from spin box)
                        self.bg1_scale = self.bg1_scale_spin.value()
                else:
                    self.bg1_scale = 1.0
                
                # Apply background subtraction to dataset 1
                i1_corrected, e1_corrected = BackgroundProcessor.subtract_background(
                    self.q1, self.i1, self.e1, self.q_bg1, self.i_bg1, self.e_bg1, self.bg1_scale)
                
                # Validate results
                warnings1 = BackgroundProcessor.validate_background_subtraction(
                    self.q1, self.i1, i1_corrected, self.bg1_scale)
                
                for warning in warnings1:
                    self.log_message(f"Dataset 1 warning: {warning}")
                
                self.log_message(f"Dataset 1: Background subtracted (scale factor: {self.bg1_scale:.4f})")
            else:
                self.log_message("Dataset 1: No background subtraction (no background file)")
            
            # Process dataset 2 if background available
            if self.q_bg2 is not None:
                # Calculate scaling factor if enabled
                if self.enable_bg_scaling.isChecked():
                    scaling_method_id = self.bg_scaling_method.checkedId()
                    
                    if scaling_method_id == 0:  # No scaling
                        self.bg2_scale = 1.0
                    elif scaling_method_id == 1:  # Auto optimize
                        self.bg2_scale = BackgroundProcessor.calculate_scaling_factor(
                            self.q2, self.i2, self.q_bg2, self.i_bg2, 'auto_optimize')
                    else:  # Manual scaling (from spin box)
                        self.bg2_scale = self.bg2_scale_spin.value()
                else:
                    self.bg2_scale = 1.0
                
                # Apply background subtraction to dataset 2
                i2_corrected, e2_corrected = BackgroundProcessor.subtract_background(
                    self.q2, self.i2, self.e2, self.q_bg2, self.i_bg2, self.e_bg2, self.bg2_scale)
                
                # Validate results
                warnings2 = BackgroundProcessor.validate_background_subtraction(
                    self.q2, self.i2, i2_corrected, self.bg2_scale)
                
                for warning in warnings2:
                    self.log_message(f"Dataset 2 warning: {warning}")
                
                self.log_message(f"Dataset 2: Background subtracted (scale factor: {self.bg2_scale:.4f})")
            else:
                self.log_message("Dataset 2: No background subtraction (no background file)")
            
            self.log_message("Background subtraction processing completed")
            
            return q1_corrected, i1_corrected, e1_corrected, q2_corrected, i2_corrected, e2_corrected
            
        except Exception as e:
            self.log_message(f"Error in background subtraction: {str(e)}")
            return self.q1, self.i1, self.e1, self.q2, self.i2, self.e2
    
    def preview_background(self):
        """Preview background subtraction results"""
        if not (self.q_bg1 is not None or self.q_bg2 is not None):
            QMessageBox.warning(self, "Missing Data", "Please load at least one background file first")
            return
        
        # Calculate current scaling factors for available backgrounds
        bg1_scale = 1.0
        bg2_scale = 1.0
        
        if self.enable_bg_scaling.isChecked():
            scaling_method_id = self.bg_scaling_method.checkedId()
            
            if scaling_method_id == 0:  # No scaling
                bg1_scale = 1.0
                bg2_scale = 1.0
            elif scaling_method_id == 1:  # Auto optimize
                if self.q_bg1 is not None:
                    bg1_scale = BackgroundProcessor.calculate_scaling_factor(
                        self.q1, self.i1, self.q_bg1, self.i_bg1, 'auto_optimize')
                    # Update spin box with calculated value
                    self.bg1_scale_spin.setValue(bg1_scale)
                if self.q_bg2 is not None:
                    bg2_scale = BackgroundProcessor.calculate_scaling_factor(
                        self.q2, self.i2, self.q_bg2, self.i_bg2, 'auto_optimize')
                    # Update spin box with calculated value
                    self.bg2_scale_spin.setValue(bg2_scale)
            else:  # Manual
                bg1_scale = self.bg1_scale_spin.value()
                bg2_scale = self.bg2_scale_spin.value()
        
        # Create dummy background data for datasets without backgrounds (for plotting)
        q_bg1_plot = self.q_bg1 if self.q_bg1 is not None else self.q1
        i_bg1_plot = self.i_bg1 if self.i_bg1 is not None else np.zeros_like(self.q1)
        e_bg1_plot = self.e_bg1 if self.e_bg1 is not None else np.zeros_like(self.q1)
        
        q_bg2_plot = self.q_bg2 if self.q_bg2 is not None else self.q2
        i_bg2_plot = self.i_bg2 if self.i_bg2 is not None else np.zeros_like(self.q2)
        e_bg2_plot = self.e_bg2 if self.e_bg2 is not None else np.zeros_like(self.q2)
        
        # If no background, set scale to 0 for visual clarity
        if self.q_bg1 is None:
            bg1_scale = 0.0
        if self.q_bg2 is None:
            bg2_scale = 0.0
        
        # Plot background comparison
        self.plot_canvas.plot_background_comparison(
            self.q1, self.i1, q_bg1_plot, i_bg1_plot, bg1_scale,
            self.q2, self.i2, q_bg2_plot, i_bg2_plot, bg2_scale,
            self.errorbar_checkbox.isChecked(), 
            self.e1, e_bg1_plot, self.e2, e_bg2_plot)
        
        # Log status
        bg1_status = f"Dataset 1: {('Background loaded, scale=' + f'{bg1_scale:.4f}') if self.q_bg1 is not None else 'No background'}"
        bg2_status = f"Dataset 2: {('Background loaded, scale=' + f'{bg2_scale:.4f}') if self.q_bg2 is not None else 'No background'}"
        self.log_message(f"Background preview: {bg1_status}, {bg2_status}")
    
    def find_auto_merge(self):
        """Find automatic merge region"""
        if self.q1 is None or self.q2 is None:
            return
            
        try:
            # Apply background subtraction if enabled
            q1_work, i1_work, e1_work, q2_work, i2_work, e2_work = self.apply_background_subtraction()
            
            self.log_message("Finding optimal merge point...")
            
            # Calculate initial merge range
            q_merge_start = self.q_overlap_start + (self.q_overlap_end - self.q_overlap_start) * 0.2
            q_merge_end = self.q_overlap_end - (self.q_overlap_end - self.q_overlap_start) * 0.2
            
            optimal_merge = SAXSProcessor.find_best_merge_point(
                q1_work, i1_work, q2_work, i2_work, q_merge_start, q_merge_end)
            
            # Use merge range around optimal point
            merge_width = (self.q_overlap_end - self.q_overlap_start) * 0.3
            q_merge_start = max(self.q_overlap_start, optimal_merge - merge_width/2)
            q_merge_end = min(self.q_overlap_end, optimal_merge + merge_width/2)
            
            # Update spinboxes
            self.q_start_spin.setValue(q_merge_start)
            self.q_end_spin.setValue(q_merge_end)
            
            self.log_message(f"Automatic merge range: {q_merge_start:.4f} to {q_merge_end:.4f}")
            
            # Preview the region
            self.preview_merge()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error finding merge region: {str(e)}")
    
    def preview_merge(self):
        """Preview the merge region"""
        if self.q1 is None or self.q2 is None:
            return
            
        q_merge_start = self.q_start_spin.value()
        q_merge_end = self.q_end_spin.value()
        
        if q_merge_start >= q_merge_end:
            QMessageBox.warning(self, "Invalid Range", "Start Q must be less than End Q")
            return
        
        # Apply background subtraction if enabled
        q1_work, i1_work, e1_work, q2_work, i2_work, e2_work = self.apply_background_subtraction()
        
        # Calculate scaling factor for preview
        scale_factor = 1.0
        if self.scaling_checkbox.isChecked():
            scale_factor = SAXSProcessor.calculate_scaling_factor(
                q1_work, i1_work, q2_work, i2_work, q_merge_start, q_merge_end, e1_work, e2_work)
        
        # Plot overlap region
        self.plot_canvas.plot_overlap_region(
            q1_work, i1_work, q2_work, i2_work, 
            self.q_overlap_start, self.q_overlap_end,
            q_merge_start, q_merge_end, scale_factor,
            self.errorbar_checkbox.isChecked(), e1_work, e2_work)
        
        self.log_message(f"Previewing merge region: {q_merge_start:.4f} to {q_merge_end:.4f}")
        if self.scaling_checkbox.isChecked():
            self.log_message(f"Scaling factor: {scale_factor:.4f}")
    
    def merge_datasets(self):
        """Perform the actual merging"""
        if self.q1 is None or self.q2 is None:
            return
            
        q_merge_start = self.q_start_spin.value()
        q_merge_end = self.q_end_spin.value()
        apply_scaling = self.scaling_checkbox.isChecked()
        
        try:
            # Apply background subtraction if enabled
            q1_work, i1_work, e1_work, q2_work, i2_work, e2_work = self.apply_background_subtraction()
            
            # Show progress
            self.progress_bar.setVisible(True)
            self.progress_bar.setRange(0, 0)  # Indeterminate progress
            
            # Create worker thread with corrected data
            self.worker = MergeWorker(q1_work, i1_work, q2_work, i2_work, 
                                     q_merge_start, q_merge_end, apply_scaling, e1_work, e2_work)
            self.worker.finished.connect(self.on_merge_finished)
            self.worker.progress.connect(self.log_message)
            self.worker.error.connect(self.on_merge_error)
            self.worker.start()
            
        except Exception as e:
            self.progress_bar.setVisible(False)
            QMessageBox.critical(self, "Error", f"Error preparing merge: {str(e)}")
    
    def on_merge_finished(self, q_merged, i_merged, e_merged, scale_factor):
        """Handle merge completion"""
        self.q_merged = q_merged
        self.i_merged = i_merged
        self.e_merged = e_merged
        self.scale_factor = scale_factor
        
        self.progress_bar.setVisible(False)
        
        # Get working data for plotting (after background subtraction)
        q1_work, i1_work, e1_work, q2_work, i2_work, e2_work = self.apply_background_subtraction()
        
        # Plot results
        self.plot_canvas.plot_final_results(
            q1_work, i1_work, q2_work, i2_work, self.q_merged, self.i_merged,
            self.q_start_spin.value(), self.q_end_spin.value(), 
            scale_factor, self.scaling_checkbox.isChecked(),
            self.errorbar_checkbox.isChecked(), e1_work, e2_work, e_merged)
        
        self.log_message(f"Merge completed! {len(q_merged)} points")
        self.log_message(f"Q range: {q_merged.min():.4f} to {q_merged.max():.4f}")
        
        if self.enable_bg_subtraction.isChecked():
            self.log_message("Background subtraction applied before merging")
        
        self.log_message("Error propagation: Properly weighted averaging in overlap region")
        
        if self.scaling_checkbox.isChecked():
            self.log_message(f"Dataset scaling factor applied: {scale_factor:.4f}")
        else:
            self.log_message("No dataset scaling applied")
        
        self.save_btn.setEnabled(True)
    
    def on_merge_error(self, error_msg):
        """Handle merge error"""
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Merge Error", f"Error during merging: {error_msg}")
    
    def save_results(self):
        """Save merged results"""
        if self.q_merged is None:
            return
        
        # Generate suggested filename
        suggested_name = self.generate_filename()
        
        # Get save location
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save merged dataset", suggested_name,
            "Data files (*.dat);;Text files (*.txt);;All files (*.*)")
        
        if file_path:
            try:
                # Save data with comprehensive metadata header
                merged_data = np.column_stack((self.q_merged, self.i_merged, self.e_merged))
                
                header_lines = [
                    f"Merged SAXS data from {os.path.basename(self.file1)} and {os.path.basename(self.file2)}",
                    f"Generated by SAXS Data Merging Tool v1.0.4",
                    f"Merge range: Q = {self.q_start_spin.value():.6f} to {self.q_end_spin.value():.6f}",
                    f"Dataset scaling applied: {'Yes' if self.scaling_checkbox.isChecked() else 'No'}",
                    f"Dataset scaling factor (dataset 2): {self.scale_factor:.6f}" if self.scaling_checkbox.isChecked() else "No dataset scaling factor applied",
                ]
                
                # Add background subtraction info
                if self.enable_bg_subtraction.isChecked():
                    header_lines.extend([
                        f"Background subtraction: Applied",
                        f"Background 1: {os.path.basename(self.bg_file1) if self.bg_file1 else 'None'}",
                        f"Background 2: {os.path.basename(self.bg_file2) if self.bg_file2 else 'None'}",
                        f"Background scaling applied: {'Yes' if self.enable_bg_scaling.isChecked() else 'No'}",
                    ])
                    if self.enable_bg_scaling.isChecked():
                        header_lines.extend([
                            f"Background 1 scaling factor: {self.bg1_scale:.6f}",
                            f"Background 2 scaling factor: {self.bg2_scale:.6f}",
                        ])
                else:
                    header_lines.append("Background subtraction: Not applied")
                
                header_lines.extend([
                    f"Error propagation: Weighted averaging in overlap region",
                    f"Total points: {len(self.q_merged)}",
                    f"Q range: {self.q_merged.min():.6f} to {self.q_merged.max():.6f}",
                    "Q(A^-1)  Intensity  Error"
                ])
                
                header = '\n'.join([f"# {line}" for line in header_lines])
                
                np.savetxt(file_path, merged_data, header=header, fmt='%.6e  %.6e  %.6e')
                
                # Save plot
                plot_path = file_path.replace('.dat', '_plot.png').replace('.txt', '_plot.png')
                self.plot_canvas.fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                
                self.log_message(f"Saved: {os.path.basename(file_path)}")
                self.log_message(f"Plot saved: {os.path.basename(plot_path)}")
                
                QMessageBox.information(self, "Save Complete", 
                                      f"Data saved to: {file_path}\nPlot saved to: {plot_path}")
                
            except Exception as e:
                QMessageBox.critical(self, "Save Error", f"Error saving file: {str(e)}")
    
    def generate_filename(self):
        """Generate smart filename suggestion"""
        def extract_meaningful_part(filename):
            base = os.path.splitext(os.path.basename(filename))[0]
            prefixes_to_remove = ['experiment', 'exp', 'data', 'sample', 'saxs', 'measurement', 'meas']
            
            for prefix in prefixes_to_remove:
                if base.lower().startswith(prefix.lower()):
                    remaining = base[len(prefix):].lstrip('_-')
                    if remaining:
                        base = remaining
                        break
            
            if len(base) > 10:
                parts = base.replace('-', '_').split('_')
                meaningful_parts = [p for p in parts if len(p) > 1]
                
                if len(meaningful_parts) >= 2:
                    base = '_'.join(meaningful_parts[:2])
                elif len(meaningful_parts) == 1:
                    base = meaningful_parts[0][:10]
                else:
                    base = base[:10]
            
            return base
        
        short1 = extract_meaningful_part(self.file1)
        short2 = extract_meaningful_part(self.file2)
        
        # Build filename suffix based on options
        suffix = "mrg"
        if self.enable_bg_subtraction.isChecked():
            suffix += "_bg"
        if self.scaling_checkbox.isChecked():
            suffix += "_s"
        
        return f"{short1}_{short2}_{suffix}.dat"

def main():
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    # Create and show main window
    window = SAXSMergeGUI()
    window.show()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()