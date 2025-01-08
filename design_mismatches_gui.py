from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,  # Added this import
    QLabel,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QSpinBox,
    QFormLayout,
    QTextEdit,
    QMessageBox,
    QApplication,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QSizePolicy,
    QAbstractScrollArea,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
import pandas as pd
import sys
import subprocess
import tempfile
import os
from io import StringIO


class PreviewFileDialog(QFileDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setOption(QFileDialog.DontUseNativeDialog, True)
        self.setViewMode(QFileDialog.Detail)

        # Create preview widget
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

        layout = self.layout()
        layout.addWidget(self.preview, 0, 3, layout.rowCount(), 1)
        self.resize(1200, 700)
        self.currentChanged.connect(self.show_preview)

    def show_preview(self, path):
        if not os.path.isfile(path):
            self.preview.clear()
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                text = "".join(f.readlines()[:50])
                self.preview.setPlainText(text)
        except Exception as e:
            self.preview.setPlainText(f"Could not preview file:\n{str(e)}")


class MismatchDesignerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.current_filtered_data = None
        self.default_params_file = os.path.join(
            os.path.dirname(__file__), "mismatch_parameters.csv"
        )
        self.generated_mismatches = None  # Store generated mismatches
        self.initUI()

    def initUI(self):
        main_layout = QVBoxLayout()

        # Back button
        back_btn = QPushButton("← Back to Main Menu")
        back_btn.clicked.connect(
            lambda: self.parent().setCurrentWidget(self.parent().widget(0))
        )
        main_layout.addWidget(back_btn)

        # Create main form layout
        form_layout = QFormLayout()

        # Add title
        design_title = QLabel("Design Mismatches")
        design_title.setStyleSheet(
            "font-size: 16px; font-weight: bold; margin: 10px 0;"
        )
        form_layout.addRow(design_title)

        # File inputs with similar styling to find_guides_gui
        self.file_input = QLineEdit()
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_file)
        file_layout = QHBoxLayout()
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(browse_btn)
        form_layout.addRow("Input TSV file:", file_layout)

        # Parameters file with default value
        self.params_input = QLineEdit()
        if os.path.exists(self.default_params_file):
            self.params_input.setText(self.default_params_file)
        params_btn = QPushButton("Browse")
        params_btn.clicked.connect(self.browse_params)
        params_layout = QHBoxLayout()
        params_layout.addWidget(self.params_input)
        params_layout.addWidget(params_btn)
        form_layout.addRow("Parameters file:", params_layout)

        # Gene input with better styling
        self.genes_input = QTextEdit()
        self.genes_input.setPlaceholderText(
            "Enter locus tags or gene names (one per line)"
        )
        self.genes_input.setMaximumHeight(100)
        form_layout.addRow("Genes:", self.genes_input)

        # Spinboxes with better styling
        self.num_variants = QSpinBox()
        self.num_variants.setRange(1, 100)
        self.num_variants.setValue(20)
        self.num_variants.setStyleSheet("padding: 5px;")
        form_layout.addRow("Number of variants:", self.num_variants)

        self.guides_per_gene = QSpinBox()
        self.guides_per_gene.setRange(1, 10)
        self.guides_per_gene.setValue(2)
        self.guides_per_gene.setStyleSheet("padding: 5px;")
        form_layout.addRow("Guides per gene:", self.guides_per_gene)

        # Generate button with matching style
        self.generate_btn = QPushButton("Generate Mismatches")
        self.generate_btn.setStyleSheet(
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
        """
        )
        self.generate_btn.clicked.connect(self.generate_mismatches)

        # Add layouts to main layout first
        main_layout.addLayout(form_layout)
        main_layout.addWidget(self.generate_btn)

        # Add save controls layout after generate button
        save_layout = QHBoxLayout()

        # Output file input
        self.output_file = QLineEdit()
        self.output_file.setPlaceholderText("Output file path")

        # Save button
        self.save_btn = QPushButton("Save")
        self.save_btn.setEnabled(False)
        self.save_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #28a745;
                color: white;
                padding: 10px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:disabled {
                background-color: #cccccc;
            }
        """
        )
        self.save_btn.clicked.connect(self.save_mismatches)

        save_layout.addWidget(self.output_file)
        save_layout.addWidget(self.save_btn)

        # Create a container widget for the save layout
        save_container = QWidget()
        save_container.setLayout(save_layout)
        main_layout.addWidget(save_container)

        # Results label
        self.results_label = QLabel("")
        self.results_label.setStyleSheet("color: #666666; margin-top: 5px;")
        main_layout.addWidget(self.results_label)

        # Add table widget
        self.table = QTableWidget()
        self.table.setStyleSheet(
            """
            QTableWidget {
                font-family: Monospace;
                font-size: 9pt;
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
        self.setWindowTitle("CRISPR Guide Mismatch Designer")

    def browse_file(self):
        dialog = PreviewFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters(["TSV Files (*.tsv)", "All Files (*)"])
        if dialog.exec_():
            chosen = dialog.selectedFiles()
            if chosen:
                input_file = chosen[0]
                self.file_input.setText(input_file)
                # Set default output filename
                default_output = input_file.replace(".tsv", "_mismatches.tsv")
                self.output_file.setText(default_output)

    def browse_params(self):
        start_dir = os.path.dirname(self.params_input.text()) or os.path.dirname(
            __file__
        )
        fname, _ = QFileDialog.getOpenFileName(
            self, "Open parameters file", start_dir, "CSV files (*.csv)"
        )
        if fname:
            self.params_input.setText(fname)

    def generate_mismatches(self):
        if not all(
            [
                self.file_input.text(),
                self.params_input.text(),
                self.genes_input.toPlainText(),
            ]
        ):
            QMessageBox.warning(self, "Error", "Please fill in all fields")
            return

        try:
            # Read input file
            data = pd.read_csv(self.file_input.text(), sep="\t")

            # Get list of genes
            genes = self.genes_input.toPlainText().strip().split("\n")

            # Filter for requested genes and sort by offset
            selected_guides = []
            num_guides = self.guides_per_gene.value()  # Get the integer value

            for gene in genes:
                gene_guides = data[
                    (
                        data["locus_tag"].str.contains(gene)
                        | data["gene"].str.contains(gene)
                    )
                    & (data["mismatches"] == 0)
                ].sort_values("offset")

                # Take the specified number of guides using integer value
                selected_guides.append(gene_guides.head(num_guides))

            # Combine all selected guides
            if selected_guides:
                selected_guides = pd.concat(selected_guides)
            else:
                raise ValueError("No matching guides found for the specified genes")

            # Use temporary file for intermediate step
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=True
            ) as tmp_file:
                selected_guides.to_csv(tmp_file.name, sep="\t", index=False)

                # Run mismatch.py
                cmd = [
                    "python",
                    "mismatch.py",
                    "--parameters_file",
                    self.params_input.text(),
                    "--num_variants",
                    str(self.num_variants.value()),
                ]

                # Capture output in memory instead of writing to file
                process = subprocess.run(
                    cmd,
                    stdin=open(tmp_file.name),
                    capture_output=True,
                    text=True,
                    check=True,
                )

            # Store the results in memory
            self.generated_mismatches = pd.read_csv(StringIO(process.stdout), sep="\t")

            # Filter out the original guides (where mismatches == 0)
            self.generated_mismatches = self.generated_mismatches[
                self.generated_mismatches["mismatches"] > 0
            ]

            # Enable save button
            self.save_btn.setEnabled(True)

            # Update table with results
            self.update_table_with_results()

            QMessageBox.information(
                self,
                "Success",
                f"Generated {len(self.generated_mismatches)} mismatched guides. Click Save to write to file.",
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def update_table_with_results(self):
        """Update the table with generated mismatches"""
        if self.generated_mismatches is None:
            return

        self.table.clear()
        self.table.setRowCount(0)
        self.table.setColumnCount(len(self.generated_mismatches.columns))
        self.table.setHorizontalHeaderLabels(self.generated_mismatches.columns)

        # Show first 1000 rows
        display_data = self.generated_mismatches.head(1000)

        for idx, row in display_data.iterrows():
            row_position = self.table.rowCount()
            self.table.insertRow(row_position)
            for col, value in enumerate(row):
                self.table.setItem(row_position, col, QTableWidgetItem(str(value)))

        # Resize columns and rows
        self.table.resizeColumnsToContents()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.resizeRowsToContents()

        # Calculate total width needed including margins
        total_width = self.table.horizontalHeader().length() + 20

        # Calculate height for minimum 10 rows plus header
        header_height = self.table.horizontalHeader().height()
        row_height = self.table.rowHeight(0) if self.table.rowCount() > 0 else 30
        min_table_height = header_height + (row_height * 10)

        # Set minimum sizes for both table and widget
        self.table.setMinimumSize(total_width, min_table_height)
        self.setMinimumWidth(total_width + 50)

        # Force layout update
        self.updateGeometry()
        if self.parent():
            self.parent().adjustSize()

        self.adjustSize()

    def save_mismatches(self):
        """Save the generated mismatches to file"""
        if self.generated_mismatches is None:
            QMessageBox.warning(self, "No Data", "No mismatches to save.")
            return

        # Get the output filename
        filename = self.output_file.text()
        if not filename:
            filename, _ = QFileDialog.getSaveFileName(
                self,
                "Save Mismatches",
                os.path.join(os.path.dirname(self.file_input.text()), "mismatches.tsv"),
                "TSV Files (*.tsv)",
            )
            if filename:
                self.output_file.setText(filename)

        if filename:
            try:
                self.generated_mismatches.to_csv(filename, sep="\t", index=False)
                QMessageBox.information(
                    self, "Success", f"Mismatches saved to:\n{filename}"
                )
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save file:\n{str(e)}")


def main():
    app = QApplication(sys.argv)
    gui = MismatchDesignerGUI()
    gui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
