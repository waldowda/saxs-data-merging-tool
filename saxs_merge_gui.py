"""
SAXS Data Merging Tool - PyQt5 GUI Version
===========================================

Version: 1.0.1
Created: 2025
Last Modified: August 2025

Description:
    A graphical tool for merging Small Angle X-ray Scattering (SAXS) datasets 
    with overlapping Q ranges. Features automatic merge region detection, 
    optional intensity scaling, and interactive visualization.

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
    - Interactive file selection
    - Automatic merge region optimization
    - Optional intensity scaling
    - Real-time preview of merge regions
    - Comprehensive data export with metadata
    - Multi-threaded processing

Citation:
    If you use this tool in your research, please cite:
    Waldow, D. (2025). SAXS Data Merging Tool (Version 1.0.0) [Computer software]. 
    Zenodo. https://doi.org/10.5281/zenodo.16734022

Changelog:
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
                           QSizePolicy)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QPixmap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class SAXSProcessor:
    """Core SAXS processing functions - same as before but as a class"""
    
    @staticmethod
    def load_saxs_data(filename, delimiter=None):
        """
        Load SAXS data from file. Assumes first column is Q, second is intensity.
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
            
            # Extract Q and intensity columns
            if data.shape[1] >= 2:
                q = data[:, 0]
                intensity = data[:, 1]
                return q, intensity
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
                return df.iloc[:, 0].values, df.iloc[:, 1].values
                
            except Exception as e2:
                raise ValueError(f"Could not load data file: {str(e2)}")  # FIXED: Removed extra )

    @staticmethod
    def find_overlap_range(q1, q2):
        q_min_overlap = max(np.min(q1), np.min(q2))
        q_max_overlap = min(np.max(q1), np.max(q2))
        
        if q_min_overlap >= q_max_overlap:
            raise ValueError("No overlap between datasets!")
        
        return q_min_overlap, q_max_overlap

    @staticmethod
    def calculate_scaling_factor(q1, i1, q2, i2, q_overlap_start, q_overlap_end):
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
    def merge_saxs_data(q1, i1, q2, i2, q_merge_start, q_merge_end, apply_scaling=True):
        scale_factor = 1.0
        i2_scaled = i2.copy()
        
        if apply_scaling:
            scale_factor = SAXSProcessor.calculate_scaling_factor(q1, i1, q2, i2, q_merge_start, q_merge_end)
            i2_scaled = i2 * scale_factor
        
        q_merged = []
        i_merged = []
        
        # Add dataset 1 up to merge start
        mask1_before = q1 < q_merge_start
        q_merged.extend(q1[mask1_before])
        i_merged.extend(i1[mask1_before])
        
        # Create overlap region with averaging
        q_overlap = np.linspace(q_merge_start, q_merge_end, 100)
        
        interp1 = interp1d(q1, i1, kind='linear', bounds_error=False, fill_value=np.nan)
        interp2 = interp1d(q2, i2_scaled, kind='linear', bounds_error=False, fill_value=np.nan)
        
        i1_overlap = interp1(q_overlap)
        i2_overlap = interp2(q_overlap)
        
        for i, q_val in enumerate(q_overlap):
            if not np.isnan(i1_overlap[i]) and not np.isnan(i2_overlap[i]):
                q_merged.append(q_val)
                i_merged.append((i1_overlap[i] + i2_overlap[i]) / 2)
            elif not np.isnan(i1_overlap[i]):
                q_merged.append(q_val)
                i_merged.append(i1_overlap[i])
            elif not np.isnan(i2_overlap[i]):
                q_merged.append(q_val)
                i_merged.append(i2_overlap[i])
        
        # Add dataset 2 after merge end
        mask2_after = q2 > q_merge_end
        q_merged.extend(q2[mask2_after])
        i_merged.extend(i2_scaled[mask2_after])
        
        q_merged = np.array(q_merged)
        i_merged = np.array(i_merged)
        
        sort_indices = np.argsort(q_merged)
        q_merged = q_merged[sort_indices]
        i_merged = i_merged[sort_indices]
        
        return q_merged, i_merged, scale_factor

class MergeWorker(QThread):
    """Worker thread for processing to avoid GUI freezing"""
    finished = pyqtSignal(object, object, float)
    progress = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, q1, i1, q2, i2, q_merge_start, q_merge_end, apply_scaling):
        super().__init__()
        self.q1, self.i1 = q1, i1
        self.q2, self.i2 = q2, i2
        self.q_merge_start = q_merge_start
        self.q_merge_end = q_merge_end
        self.apply_scaling = apply_scaling
    
    def run(self):
        try:
            self.progress.emit("Merging datasets...")
            q_merged, i_merged, scale_factor = SAXSProcessor.merge_saxs_data(
                self.q1, self.i1, self.q2, self.i2, 
                self.q_merge_start, self.q_merge_end, self.apply_scaling)
            self.finished.emit(q_merged, i_merged, scale_factor)
        except Exception as e:
            self.error.emit(str(e))

class PlotCanvas(FigureCanvas):
    """Matplotlib canvas for displaying plots in Qt"""
    
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(12, 8))
        super().__init__(self.fig)
        self.setParent(parent)
        
    def plot_overlap_region(self, q1, i1, q2, i2, q_overlap_start, q_overlap_end, 
                           q_merge_start=None, q_merge_end=None, scale_factor=1.0):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        
        # Plot full datasets (FIXED: Show complete Q range)
        ax.loglog(q1, i1, 'b-', label='Dataset 1', linewidth=2, alpha=0.7)
        ax.loglog(q2, i2, 'r-', label='Dataset 2 (original)', alpha=0.5)
        ax.loglog(q2, i2 * scale_factor, 'r--', label='Dataset 2 (scaled)', linewidth=2)
        
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
                          q_merge_start, q_merge_end, scale_factor, apply_scaling):
        self.fig.clear()
        
        # Original datasets
        ax1 = self.fig.add_subplot(2, 1, 1)
        ax1.loglog(q1, i1, 'b-', label='Dataset 1', alpha=0.7)
        ax1.loglog(q2, i2, 'r-', label='Dataset 2', alpha=0.7)
        if apply_scaling:
            ax1.loglog(q2, i2 * scale_factor, 'r--', label='Dataset 2 (scaled)', alpha=0.7)
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
        ax2.loglog(q_merged, i_merged, 'k-', linewidth=2, label='Merged dataset')
        ax2.loglog(q1, i1, 'b-', alpha=0.3, label='Dataset 1')
        if apply_scaling:
            ax2.loglog(q2, i2 * scale_factor, 'r-', alpha=0.3, label='Dataset 2 (scaled)')
        else:
            ax2.loglog(q2, i2, 'r-', alpha=0.3, label='Dataset 2')
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
        self.q1 = self.i1 = self.q2 = self.i2 = None
        self.q_overlap_start = self.q_overlap_end = None
        self.q_merged = self.i_merged = None
        self.scale_factor = 1.0
        self.file1 = self.file2 = ""
        
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('SAXS Data Merging Tool')
        self.setGeometry(100, 100, 1200, 900)
        
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
        version_label = QLabel('Version 1.0.0')
        version_label.setFont(QFont('Arial', 10))
        version_label.setAlignment(Qt.AlignCenter)
        version_label.setStyleSheet("color: gray;")
        title_layout.addWidget(version_label)
        
        # Contact/License info space - customize as needed
        info_label = QLabel(
            'Author: Dean Waldow | Email: waldowda@plu.edu\n'
            'Institution: Pacific Lutheran University | License: MIT\n'
            'DOI: 10.5281/zenodo.16734022 | GitHub: github.com/waldowda/saxs-data-merging-tool'
        )
        info_label.setFont(QFont('Arial', 9))
        info_label.setAlignment(Qt.AlignCenter)
        info_label.setStyleSheet("color: #666666; padding: 5px;")
        info_label.setWordWrap(True)
        title_layout.addWidget(info_label)
        
        # Add some spacing
        title_layout.addSpacing(10)
        
        layout.addWidget(title_widget)
        
        # File selection group
        file_group = QGroupBox("File Selection")
        file_layout = QGridLayout(file_group)
        
        self.file1_label = QLabel("No file selected")
        self.file2_label = QLabel("No file selected")
        
        file_layout.addWidget(QLabel("Dataset 1 (lower Q data):"), 0, 0)
        file_layout.addWidget(self.file1_label, 0, 1)
        file_layout.addWidget(QPushButton("Browse...", clicked=self.select_file1), 0, 2)
        
        file_layout.addWidget(QLabel("Dataset 2 (higher Q data):"), 1, 0)
        file_layout.addWidget(self.file2_label, 1, 1)
        file_layout.addWidget(QPushButton("Browse...", clicked=self.select_file2), 1, 2)
        
        layout.addWidget(file_group)
        
        # Options group
        options_group = QGroupBox("Merge Options")
        options_layout = QVBoxLayout(options_group)
        
        self.scaling_checkbox = QCheckBox("Apply intensity scaling to match datasets")
        self.scaling_checkbox.setChecked(True)
        options_layout.addWidget(self.scaling_checkbox)
        
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
        
        self.new_analysis_btn = QPushButton("New Analysis", clicked=self.reset_analysis)
        self.merge_btn = QPushButton("Merge Datasets", clicked=self.merge_datasets)
        self.merge_btn.setEnabled(False)
        
        self.save_btn = QPushButton("Save Results", clicked=self.save_results)
        self.save_btn.setEnabled(False)
        
        button_layout.addWidget(self.new_analysis_btn)
        button_layout.addWidget(self.merge_btn)
        button_layout.addWidget(self.save_btn)
        button_layout.addStretch()
        
        layout.addLayout(button_layout)
        
        # Status and progress
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
    def log_message(self, message):
        """Add message to status text"""
        self.status_text.append(message)
        QApplication.processEvents()
    
    def reset_analysis(self):
        """Reset all analysis data and UI state for new analysis"""
        # Stop any running worker
        if hasattr(self, 'worker') and self.worker.isRunning():
            self.worker.terminate()
            self.worker.wait()
        
        # Clear data
        self.q1 = self.i1 = self.q2 = self.i2 = None
        self.q_overlap_start = self.q_overlap_end = None
        self.q_merged = self.i_merged = None
        self.scale_factor = 1.0
        
        # Reset file paths and labels
        self.file1 = self.file2 = ""
        self.file1_label.setText("No file selected")
        self.file2_label.setText("No file selected")
        
        # Reset UI state
        self.q_start_spin.setRange(0.0, 10.0)
        self.q_end_spin.setRange(0.0, 10.0)
        self.q_start_spin.setValue(0.0)
        self.q_end_spin.setValue(0.0)
        
        # Reset scaling checkbox to default
        self.scaling_checkbox.setChecked(True)
        
        # Disable controls
        self.auto_merge_btn.setEnabled(False)
        self.preview_btn.setEnabled(False)
        self.merge_btn.setEnabled(False)
        self.save_btn.setEnabled(False)
        
        # Clear plot
        self.plot_canvas.fig.clear()
        self.plot_canvas.draw()
        
        # Hide progress bar
        self.progress_bar.setVisible(False)
        
        # Clear status log
        self.status_text.clear()
        self.log_message("--- New Analysis Started ---")
        
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
    
    def load_and_analyze(self):
        """Load data and analyze overlap if both files selected"""
        if not (self.file1 and self.file2):
            return
            
        try:
            # Load data and capture Xenocs header message (FIXED: Log to GUI)
            self.log_message("Loading SAXS datasets...")
            self.q1, self.i1 = SAXSProcessor.load_saxs_data(self.file1)
            
            # Check if Xenocs header was detected and log to GUI
            if self.file1:
                with open(self.file1, 'r') as f:
                    first_line = f.readline().strip()
                if first_line.startswith('###########'):
                    self.log_message("Detected Xenocs header in dataset 1 - automatically skipped")
            
            self.q2, self.i2 = SAXSProcessor.load_saxs_data(self.file2)
            
            # Check if Xenocs header was detected and log to GUI
            if self.file2:
                with open(self.file2, 'r') as f:
                    first_line = f.readline().strip()
                if first_line.startswith('###########'):
                    self.log_message("Detected Xenocs header in dataset 2 - automatically skipped")
            
            self.log_message(f"Dataset 1: Q range {self.q1.min():.4f} to {self.q1.max():.4f} ({len(self.q1)} points)")
            self.log_message(f"Dataset 2: Q range {self.q2.min():.4f} to {self.q2.max():.4f} ({len(self.q2)} points)")
            
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
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading data: {str(e)}")
    
    def find_auto_merge(self):
        """Find automatic merge region"""
        if self.q1 is None or self.q2 is None:
            return
            
        try:
            self.log_message("Finding optimal merge point...")
            
            # Calculate initial merge range
            q_merge_start = self.q_overlap_start + (self.q_overlap_end - self.q_overlap_start) * 0.2
            q_merge_end = self.q_overlap_end - (self.q_overlap_end - self.q_overlap_start) * 0.2
            
            optimal_merge = SAXSProcessor.find_best_merge_point(
                self.q1, self.i1, self.q2, self.i2, q_merge_start, q_merge_end)
            
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
        
        # Calculate scaling factor for preview
        scale_factor = 1.0
        if self.scaling_checkbox.isChecked():
            scale_factor = SAXSProcessor.calculate_scaling_factor(
                self.q1, self.i1, self.q2, self.i2, q_merge_start, q_merge_end)
        
        # Plot overlap region
        self.plot_canvas.plot_overlap_region(
            self.q1, self.i1, self.q2, self.i2, 
            self.q_overlap_start, self.q_overlap_end,
            q_merge_start, q_merge_end, scale_factor)
        
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
        
        # Show progress
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Create worker thread
        self.worker = MergeWorker(self.q1, self.i1, self.q2, self.i2, 
                                 q_merge_start, q_merge_end, apply_scaling)
        self.worker.finished.connect(self.on_merge_finished)
        self.worker.progress.connect(self.log_message)
        self.worker.error.connect(self.on_merge_error)
        self.worker.start()
    
    def on_merge_finished(self, q_merged, i_merged, scale_factor):
        """Handle merge completion"""
        self.q_merged = q_merged
        self.i_merged = i_merged
        self.scale_factor = scale_factor
        
        self.progress_bar.setVisible(False)
        
        # Plot results
        self.plot_canvas.plot_final_results(
            self.q1, self.i1, self.q2, self.i2, self.q_merged, self.i_merged,
            self.q_start_spin.value(), self.q_end_spin.value(), 
            scale_factor, self.scaling_checkbox.isChecked())
        
        self.log_message(f"Merge completed! {len(q_merged)} points")
        self.log_message(f"Q range: {q_merged.min():.4f} to {q_merged.max():.4f}")
        
        if self.scaling_checkbox.isChecked():
            self.log_message(f"Scaling factor applied: {scale_factor:.4f}")
        else:
            self.log_message("No scaling applied")
        
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
                # Save data with metadata header
                merged_data = np.column_stack((self.q_merged, self.i_merged))
                
                header_lines = [
                    f"Merged SAXS data from {os.path.basename(self.file1)} and {os.path.basename(self.file2)}",
                    f"Merge range: Q = {self.q_start_spin.value():.6f} to {self.q_end_spin.value():.6f}",
                    f"Scaling applied: {'Yes' if self.scaling_checkbox.isChecked() else 'No'}",
                    f"Scaling factor (dataset 2): {self.scale_factor:.6f}" if self.scaling_checkbox.isChecked() else "No scaling factor applied",
                    f"Total points: {len(self.q_merged)}",
                    f"Q range: {self.q_merged.min():.6f} to {self.q_merged.max():.6f}",
                    "Q(A^-1)  Intensity"
                ]
                header = '\n'.join([f"# {line}" for line in header_lines])
                
                np.savetxt(file_path, merged_data, header=header, fmt='%.6e  %.6e')
                
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
        
        if self.scaling_checkbox.isChecked():
            return f"{short1}_{short2}_mrg_s.dat"
        else:
            return f"{short1}_{short2}_mrg.dat"

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