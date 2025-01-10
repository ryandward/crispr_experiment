import sys
import os
import argparse
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QComboBox,
    QMessageBox,
    QProgressDialog,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QFormLayout,
    QHBoxLayout,
    QSizePolicy,
    QCheckBox,
    QSpinBox,
    QAbstractScrollArea,
    QGridLayout,
    QTextEdit,
)
from PyQt5.QtCore import Qt, QRegExp, QProcess
from PyQt5.QtGui import QRegExpValidator, QColor, QPixmap
import csv
import pandas as pd
from io import StringIO


class PreviewFileDialog(QFileDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setOption(QFileDialog.DontUseNativeDialog, True)
        self.setViewMode(QFileDialog.Detail)

        # Create preview widget with better styling
        self.preview = QTextEdit(self)
        self.preview.setReadOnly(True)
        self.preview.setMinimumWidth(400)
        self.preview.setStyleSheet(
            """
            QTextEdit {
                font-family: Monospace;
                font-size: 9pt;
                background-color: #ffffff;
                padding: 5px;
            }
        """
        )

        # Get the layout and add preview to the right side
        layout = self.layout()
        layout.addWidget(self.preview, 0, 3, layout.rowCount(), 1)

        # Set better initial size
        self.resize(1200, 700)

        self.currentChanged.connect(self.show_preview)

    def show_preview(self, path):
        if not os.path.isfile(path):
            self.preview.clear()
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                # Read first 50 lines for better preview
                text = "".join(f.readlines()[:50])
                self.preview.setPlainText(text)
        except Exception as e:
            self.preview.setPlainText(f"Could not preview file:\n{str(e)}")


class FindGuidesGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.process = None
        self.all_data = None  # Will be a pandas DataFrame
        self.stdout_buffer = StringIO()  # Buffer to collect all stdout
        self.headers = []  # Store headers
        self.filter_widget = None  # Add this
        self.initUI()

    def initUI(self):
        main_layout = QVBoxLayout()

        back_btn = QPushButton("← Back to Main Menu")
        back_btn.clicked.connect(
            lambda: self.parent().setCurrentWidget(
                self.parent().widget(0)  # Assuming welcome page is at index 0
            )
        )
        main_layout.addWidget(back_btn)

        # Create a widget for the "Find all Guides" section
        find_guides_widget = QWidget()
        find_guides_layout = QFormLayout()

        # Add Find Guides Section Title
        find_title = QLabel("Find all Guides")
        find_title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px 0;")
        find_guides_layout.addRow(find_title)

        # Group all the input controls
        self.genome_file_edit = QLineEdit()
        browse_genome_btn = QPushButton("Browse")
        browse_genome_btn.clicked.connect(self.browse_genome_file)
        genome_layout = QHBoxLayout()
        genome_layout.addWidget(self.genome_file_edit)
        genome_layout.addWidget(browse_genome_btn)
        find_guides_layout.addRow(QLabel("Genome File"), genome_layout)

        self.pam_edit = QLineEdit("NGG")
        find_guides_layout.addRow(QLabel("PAM Sequence"), self.pam_edit)

        self.barcode_edit = QLineEdit("20")
        find_guides_layout.addRow(QLabel("Barcode Length"), self.barcode_edit)

        self.mismatches_combo = QComboBox()
        self.mismatches_combo.addItems(["0", "1", "2", "3"])  # Add "3" to options
        self.mismatches_combo.setCurrentText("3")  # Set default to "3"
        find_guides_layout.addRow(QLabel("Off-target Tolerance"), self.mismatches_combo)

        self.pam_direction_combo = QComboBox()
        self.pam_direction_combo.addItems(["upstream", "downstream"])
        self.pam_direction_combo.setCurrentText("downstream")  # Pre-select downstream
        find_guides_layout.addRow(QLabel("PAM Direction"), self.pam_direction_combo)

        self.regulatory_edit = QLineEdit("0")
        find_guides_layout.addRow(QLabel("Regulatory Size"), self.regulatory_edit)

        # Add the run button and results label to the find guides section
        self.run_button = QPushButton(
            "Find all possible guides"
        )  # Make button accessible as instance variable
        self.run_button.setStyleSheet(
            """
            QPushButton {
                background-color: #007bff;
                color: white;
                padding: 10px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #0056b3;
            }
            QPushButton:disabled {
                background-color: #cccccc;
                color: #666666;
            }
        """
        )
        self.run_button.clicked.connect(self.run_script)

        self.results_label = QLabel("")  # Add label for showing number of guides
        self.results_label.setStyleSheet("color: #666666; margin-top: 5px;")

        find_guides_layout.addRow(self.run_button)
        find_guides_layout.addRow(self.results_label)

        find_guides_widget.setLayout(find_guides_layout)
        main_layout.addWidget(find_guides_widget)

        # Table setup
        self.table = QTableWidget()
        self.table.setStyleSheet(
            """
            QTableWidget {
                font-family: Monospace;
                font-size: 9pt;
                color: #000000;

            }
            QHeaderView::section {
                font-family: Monospace;
                font-size: 9pt;
                font-weight: bold;
            }
        """
        )
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.table.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        main_layout.addWidget(self.table)

        self.setLayout(main_layout)

    def browse_genome_file(self):
        dialog = PreviewFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters(["GenBank Files (*.gb *.gbk)", "All Files (*)"])
        if dialog.exec_():
            chosen = dialog.selectedFiles()
            if chosen:
                self.genome_file_edit.setText(chosen[0])

    def run_script(self):
        genome_file = self.genome_file_edit.text().strip()
        pam = self.pam_edit.text().strip()
        barcode_length = self.barcode_edit.text().strip()
        if not genome_file or not pam or not barcode_length:
            QMessageBox.warning(self, "Input Error", "Please fill in all fields.")
            return

        self.table.setRowCount(0)
        script_path = os.path.join(os.path.dirname(__file__), "find_guides.py")
        mismatches = self.mismatches_combo.currentText().strip()
        pam_direction = self.pam_direction_combo.currentText()
        regulatory = self.regulatory_edit.text().strip()

        args = [
            "--genome-file",
            genome_file,
            "--pam",
            pam,
            "--barcode-length",
            barcode_length,
            "--mismatches",
            mismatches,
            "--pam-direction",
            pam_direction,
            "--regulatory",
            regulatory,
        ]

        self.progress = QProgressDialog("Running find_guides...", "Cancel", 0, 0, self)
        self.progress.setWindowModality(Qt.WindowModal)
        self.progress.canceled.connect(self.cancel_process)
        self.progress.show()

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.SeparateChannels)
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.process_finished)
        self.process.start(sys.executable, [script_path] + args)

    def handle_stdout(self):
        output = self.process.readAllStandardOutput().data().decode("utf-8")
        self.stdout_buffer.write(output)  # Collect all output

    def process_finished(self):
        self.progress.close()

        # Now parse the complete output as a dataframe
        self.stdout_buffer.seek(0)
        self.all_data = pd.read_csv(self.stdout_buffer, sep="\t")
        self.stdout_buffer.close()

        # Sort the data by locus_tag and offset
        self.all_data = self.all_data.sort_values(["locus_tag", "offset"])

        # Get the max offset value from the data
        max_offset = self.all_data["offset"].max()

        # Set up table headers
        self.headers = list(self.all_data.columns)
        self.table.setColumnCount(len(self.headers))
        self.table.setHorizontalHeaderLabels(self.headers)

        # Create and show filter controls
        if not self.filter_widget:
            barcode_length = int(self.barcode_edit.text().strip())
            self.filter_widget = QWidget()
            filter_layout = QFormLayout()

            # Add Filter Section Title
            filter_title = QLabel("Filter Guides")
            filter_title.setStyleSheet(
                "font-size: 16px; font-weight: bold; margin: 10px 0;"
            )
            filter_layout.addRow(filter_title)

            # Create all filter controls but DON'T connect signals yet
            self.orientation_filter_combo = QComboBox()
            self.orientation_filter_combo.addItems(
                [
                    "Any",
                    "Same (Spacer Orientation == Target Orientation)",
                    "Opposite (Spacer Orientation != Target Orientation)",
                ]
            )

            self.offtarget_filter = QCheckBox("No off-targets (Limit to 1 site)")
            self.no_overlapping_guides_filter = QCheckBox("No overlapping guides")
            self.no_overlapping_genes_filter = QCheckBox(
                "No overlapping genes"
            )  # Add new checkbox

            # Create horizontal layout for checkboxes
            checkbox_layout = QHBoxLayout()
            checkbox_layout.addWidget(self.offtarget_filter)
            checkbox_layout.addWidget(self.no_overlapping_guides_filter)
            checkbox_layout.addWidget(self.no_overlapping_genes_filter)  # Add to layout
            filter_layout.addRow(checkbox_layout)

            self.no_min_overlap = QCheckBox("No minimum")
            self.no_min_overlap.stateChanged.connect(self.toggle_min_overlap)
            self.min_overlap = QSpinBox()
            self.min_overlap.setRange(0, barcode_length)
            self.min_overlap.setValue(barcode_length)  # Set default to barcode length

            self.no_max_overlap = QCheckBox("No maximum")
            self.no_max_overlap.stateChanged.connect(self.toggle_max_overlap)
            self.max_overlap = QSpinBox()
            self.max_overlap.setRange(0, barcode_length)
            self.max_overlap.setValue(barcode_length)  # Set default to barcode length

            overlap_layout = QHBoxLayout()
            overlap_layout.addWidget(QLabel("Min:"))
            overlap_layout.addWidget(self.min_overlap)
            overlap_layout.addWidget(self.no_min_overlap)
            overlap_layout.addWidget(QLabel("Max:"))
            overlap_layout.addWidget(self.max_overlap)
            overlap_layout.addWidget(self.no_max_overlap)

            # Add Offset range filter with No min/max options
            self.no_min_offset = QCheckBox("No minimum")
            self.no_min_offset.stateChanged.connect(self.toggle_min_offset)
            self.min_offset = QSpinBox()
            self.min_offset.setRange(-999999, 999999)

            self.no_max_offset = QCheckBox("No maximum")
            self.no_max_offset.stateChanged.connect(self.toggle_max_offset)
            self.max_offset = QSpinBox()
            self.max_offset.setRange(-999999, 999999)
            self.max_offset.setValue(max_offset)

            offset_layout = QHBoxLayout()
            offset_layout.addWidget(QLabel("Min:"))
            offset_layout.addWidget(self.min_offset)
            offset_layout.addWidget(self.no_min_offset)
            offset_layout.addWidget(QLabel("Max:"))
            offset_layout.addWidget(self.max_offset)
            offset_layout.addWidget(self.no_max_offset)

            # Add Apply Filters button
            self.apply_filters_btn = QPushButton("Apply Filters")
            self.apply_filters_btn.setStyleSheet(
                """
                QPushButton {
                    background-color: #28a745;
                    color: white;
                    padding: 8px;
                    border-radius: 4px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #218838;
                }
            """
            )
            self.apply_filters_btn.clicked.connect(self.apply_filters)

            self.save_button = QPushButton("Save As...")
            self.save_button.setStyleSheet(
                """
                QPushButton {
                    background-color: #17a2b8;
                    color: white;
                    padding: 8px;
                    border-radius: 4px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #138496;
                }
                """
            )
            self.save_button.clicked.connect(self.save_filtered_data)

            # Add all widgets to layout
            filter_layout.addRow(
                QLabel("Guide/Gene Orientation"), self.orientation_filter_combo
            )
            filter_layout.addRow("Gene Body Overlap:", overlap_layout)
            filter_layout.addRow("Gene Body Offset:", offset_layout)

            # Add guides per gene filter
            self.guides_per_gene = QSpinBox()
            self.guides_per_gene.setRange(0, 100)
            self.guides_per_gene.setValue(0)  # 0 means no limit
            self.guides_per_gene.setSpecialValueText(
                "All"
            )  # Show "All" when value is 0
            filter_layout.addRow("Guides per gene:", self.guides_per_gene)

            filter_layout.addRow(self.apply_filters_btn)
            filter_layout.addRow(self.save_button)

            # Add filter results label
            self.filter_results_label = QLabel("")
            self.filter_results_label.setStyleSheet("color: #666666; margin-top: 5px;")
            filter_layout.addRow(self.filter_results_label)

            self.filter_widget.setLayout(filter_layout)
            self.layout().insertWidget(
                self.layout().indexOf(self.table), self.filter_widget
            )

        # Show initial data and update counts
        self.display_all_data()
        self.run_button.setEnabled(False)  # Disable the run button
        self.results_label.setText(f"Found {len(self.all_data):,} guides")

        QMessageBox.information(
            self,
            "Process Finished",
            f"Loaded {len(self.all_data):,} guides. Use filters and click Apply to explore results.",
        )

    def toggle_min_offset(self, state):
        self.min_offset.setEnabled(not state)

    def toggle_max_offset(self, state):
        self.max_offset.setEnabled(not state)

    def toggle_min_overlap(self, state):
        self.min_overlap.setEnabled(not state)

    def toggle_max_overlap(self, state):
        self.max_overlap.setEnabled(not state)

    def display_all_data(self):
        """Display first 1000 rows of sorted data"""
        self.table.setRowCount(0)

        # Group by locus_tag and take first 1000 rows
        displayed_data = self.all_data.head(1000)

        current_locus = None
        use_alternate_color = False

        for _, row in displayed_data.iterrows():
            table_row = self.table.rowCount()
            self.table.insertRow(table_row)

            # Change background color when locus_tag changes
            if current_locus != row["locus_tag"]:
                current_locus = row["locus_tag"]
                use_alternate_color = not use_alternate_color

            background_color = "#e6f3ff" if use_alternate_color else "#ffffff"

            for col_idx, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                item.setBackground(QColor(background_color))
                self.table.setItem(table_row, col_idx, item)

        self.table.resizeColumnsToContents()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.resizeRowsToContents()

        # Calculate total width needed including margins with extra padding
        total_width = self.table.horizontalHeader().length() + 20

        # Calculate height for minimum 10 rows plus header
        header_height = self.table.horizontalHeader().height()
        row_height = self.table.rowHeight(0)  # Get height of first row
        min_table_height = header_height + (row_height * 10)

        # Set minimum sizes for both table and widget
        self.table.setMinimumSize(total_width, min_table_height)
        self.setMinimumWidth(total_width + 50)  # Additional padding for widget

        # Force layout update
        self.updateGeometry()
        if self.parent():
            self.parent().adjustSize()

        self.adjustSize()

    def apply_filters(self):
        if self.all_data is None:
            return

        # Create progress dialog for filtering
        progress = QProgressDialog("Applying filters...", "Cancel", 0, 100, self)
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()

        self.table.setRowCount(0)
        filtered_data = self.all_data.copy()
        progress.setValue(10)

        # Apply Orientation Filter
        filter_choice = self.orientation_filter_combo.currentText()
        if filter_choice.startswith("Same"):
            filtered_data = filtered_data[
                filtered_data["sp_dir"] == filtered_data["tar_dir"]
            ]
        elif filter_choice.startswith("Opposite"):
            filtered_data = filtered_data[
                filtered_data["sp_dir"] != filtered_data["tar_dir"]
            ]
        progress.setValue(20)

        # Apply Off-target Filter
        if self.offtarget_filter.isChecked():
            filtered_data = filtered_data[filtered_data["sites"] == 1]
        progress.setValue(30)

        # Apply no overlapping genes filter
        if self.no_overlapping_genes_filter.isChecked():
            filtered_data = filtered_data[filtered_data["genes"] <= 1]
        progress.setValue(35)

        # Apply overlap filter
        progress.setLabelText("Applying overlap filters...")
        min_overlap = self.min_overlap.value()
        max_overlap = self.max_overlap.value()
        filtered_data = filtered_data[
            (filtered_data["overlap"] >= min_overlap)
            & (filtered_data["overlap"] <= max_overlap)
        ]
        progress.setValue(40)

        # Apply offset filter
        progress.setLabelText("Applying offset filters...")
        if not self.no_min_offset.isChecked():
            filtered_data = filtered_data[
                filtered_data["offset"] >= self.min_offset.value()
            ]
        if not self.no_max_offset.isChecked():
            filtered_data = filtered_data[
                filtered_data["offset"] <= self.max_offset.value()
            ]
        progress.setValue(50)

        # Apply no overlapping guides filter - now as the last step
        if self.no_overlapping_guides_filter.isChecked():
            progress.setLabelText("Processing non-overlapping guides...")
            non_overlapping = []
            total_groups = len(filtered_data.groupby("locus_tag"))
            for idx, (_, group) in enumerate(filtered_data.groupby("locus_tag")):
                sorted_guides = group.sort_values("offset")
                current_end = float("-inf")

                for _, guide in sorted_guides.iterrows():
                    guide_length = int(self.barcode_edit.text().strip())
                    guide_end = guide["offset"] + guide_length
                    if guide["offset"] >= current_end:
                        non_overlapping.append(guide)
                        current_end = guide_end

                progress.setValue(50 + (40 * idx // total_groups))
                if progress.wasCanceled():
                    return

            filtered_data = pd.DataFrame(non_overlapping)
        progress.setValue(90)

        # Apply guides per gene filter as the final step
        guides_limit = self.guides_per_gene.value()
        if guides_limit > 0:
            progress.setLabelText("Applying guides per gene limit...")
            filtered_groups = []
            total_groups = len(filtered_data.groupby("locus_tag"))

            for idx, (_, group) in enumerate(filtered_data.groupby("locus_tag")):
                filtered_groups.append(group.head(guides_limit))
                progress.setValue(90 + (8 * idx // total_groups))
                if progress.wasCanceled():
                    return

            filtered_data = pd.concat(filtered_groups)

        # Sort and store filtered data
        progress.setLabelText("Finalizing results...")
        filtered_data = filtered_data.sort_values(["locus_tag", "offset"])
        self.current_filtered_data = filtered_data

        # Update table with first 1000 rows
        self.update_table_with_data(filtered_data.head(1000))
        progress.setValue(100)

        # Update filter results count
        self.filter_results_label.setText(
            f"Showing {len(filtered_data):,} guides after filtering"
        )

        # Show completion message
        QMessageBox.information(
            self,
            "Filter Applied",
            f"Showing first 1000 of {len(filtered_data):,} matching guides",
        )

    def update_table_with_data(self, data):
        """Helper method to update table with filtered data"""
        current_locus = None
        use_alternate_color = False

        for _, row in data.iterrows():
            table_row = self.table.rowCount()
            self.table.insertRow(table_row)

            if current_locus != row["locus_tag"]:
                current_locus = row["locus_tag"]
                use_alternate_color = not use_alternate_color

            bg_color = "#e6f3ff" if use_alternate_color else "#ffffff"

            for col_idx, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                item.setBackground(QColor(bg_color))
                self.table.setItem(table_row, col_idx, item)

    def save_filtered_data(self):
        if (
            not hasattr(self, "current_filtered_data")
            or self.current_filtered_data is None
        ):
            QMessageBox.warning(self, "No Data", "No filtered data to save.")
            return

        fname, _ = QFileDialog.getSaveFileName(
            self,
            "Save Filtered Data",
            os.path.join(os.path.dirname(__file__), "filtered_data.tsv"),
            "TSV Files (*.tsv)",
        )
        if not fname:
            return

        self.current_filtered_data.to_csv(fname, sep="\t", index=False)
        QMessageBox.information(self, "Data Saved", f"Filtered data saved to:\n{fname}")

    def handle_stderr(self):
        stderr = self.process.readAllStandardError().data().decode()
        if stderr:
            print(stderr)

    def cancel_process(self):
        if self.process and self.process.state() != QProcess.NotRunning:
            self.process.kill()


def main():
    from PyQt5.QtWidgets import QApplication

    app = QApplication(sys.argv)
    gui = FindGuidesGUI()
    gui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
