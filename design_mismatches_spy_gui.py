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
    QProgressDialog,  # Added this import
    QComboBox,  # Added this import
)
from PyQt5.QtCore import Qt, QProcess  # Added QProcess import
from PyQt5.QtGui import QColor
import pandas as pd
import sys
import subprocess
import tempfile
import os
from io import StringIO
from Bio.Seq import Seq  # Add this import
from rich.table import Table
from rich.console import Console
from rich.box import DOUBLE_EDGE  # Add this import


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
        self.generated_mismatches = None
        self.off_target_results = None  # New attribute to store off-target results
        self.process = None  # Add QProcess instance variable
        self.console = Console(file=sys.stderr)  # Add this line
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

        # Remove the old layout for genes input and "Load from File" button

        # Create QLineEdit and "Browse" button for genes file path
        self.genes_file_input = QLineEdit()
        genes_browse_btn = QPushButton("Browse")
        genes_browse_btn.clicked.connect(self.load_genes_with_preview)

        genes_path_layout = QHBoxLayout()
        genes_path_layout.addWidget(self.genes_file_input)
        genes_path_layout.addWidget(genes_browse_btn)
        form_layout.addRow("Genes or locus tags:", genes_path_layout)

        # Add the existing QTextEdit on a separate row for manual editing
        form_layout.addRow("", self.genes_input)

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

        # Off-target Settings Section
        self.off_target_settings_box = QWidget()
        self.off_target_settings_box.setEnabled(False)
        off_target_layout = QFormLayout(self.off_target_settings_box)

        # Add disabled state styling
        self.off_target_settings_box.setStyleSheet(
            """
            QWidget:disabled {
                color: #999999;
            }
            QLineEdit:disabled, QComboBox:disabled {
                background-color: #f0f0f0;
                color: #999999;
            }
        """
        )

        off_target_settings_title = QLabel(
            "Off-target Settings (Available after generating mismatches)"
        )
        off_target_settings_title.setStyleSheet(
            """
            font-size: 14px;
            font-weight: bold;
            margin: 10px 0;
            color: #999999;
            """
        )
        off_target_layout.addRow(off_target_settings_title)

        self.pam_edit = QLineEdit("NGG")
        self.pam_edit.setPlaceholderText("Available after generating mismatches")
        off_target_layout.addRow("PAM Sequence:", self.pam_edit)

        self.mismatches_combo = QComboBox()
        self.mismatches_combo.addItems(["0", "1", "2", "3"])  # Add "3" to options
        self.mismatches_combo.setCurrentText("3")  # Set default to "3"
        off_target_layout.addRow("Off-target Tolerance:", self.mismatches_combo)

        self.pam_direction_combo = QComboBox()
        self.pam_direction_combo.addItems(["upstream", "downstream"])
        self.pam_direction_combo.setCurrentText("downstream")
        off_target_layout.addRow("PAM Direction:", self.pam_direction_combo)

        main_layout.addWidget(self.off_target_settings_box)

        # Genome file section
        genome_layout = QHBoxLayout()
        self.genome_input = QLineEdit()
        self.genome_input.setPlaceholderText(
            "Select genome file (enabled after generating mismatches)"
        )
        self.genome_input.setEnabled(False)

        # Make browse_genome_btn a class attribute
        self.browse_genome_btn = QPushButton("Browse")
        self.browse_genome_btn.clicked.connect(self.browse_genome_file)
        self.browse_genome_btn.setEnabled(False)

        genome_layout.addWidget(self.genome_input)
        genome_layout.addWidget(self.browse_genome_btn)

        # Style disabled genome components
        self.genome_input.setStyleSheet(
            """
            QLineEdit:disabled {
                background-color: #f0f0f0;
                color: #999999;
            }
        """
        )

        self.browse_genome_btn.setStyleSheet(
            """
            QPushButton:disabled {
                background-color: #f0f0f0;
                color: #999999;
            }
        """
        )

        # Check Off-targets button
        self.check_offtargets_btn = QPushButton("Remove off-targets")
        self.check_offtargets_btn.setEnabled(False)
        self.check_offtargets_btn.clicked.connect(
            self.check_offtargets
        )  # Add this line
        self.check_offtargets_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #17a2b8;
                color: white;
                padding: 10px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #138496;
            }
            QPushButton:disabled {
                background-color: #f0f0f0;
                color: #999999;
                border: 1px solid #cccccc;
            }
        """
        )

        genome_layout.addWidget(self.check_offtargets_btn)

        # Create container for genome selection
        genome_container = QWidget()
        genome_container.setLayout(genome_layout)
        main_layout.addWidget(genome_container)

        # Move save controls here (after genome container)
        save_layout = QHBoxLayout()
        self.output_file = QLineEdit()
        self.output_file.setPlaceholderText("Output file path")
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

            # Calculate total operations for progress bar
            total_guides = len(genes) * self.guides_per_gene.value()
            total_operations = total_guides * self.num_variants.value()

            # Create and configure progress dialog
            progress = QProgressDialog(
                "Generating mismatches...", "Cancel", 0, total_operations, self
            )
            progress.setWindowModality(Qt.WindowModal)
            progress.setMinimumDuration(0)
            progress.setValue(0)

            # Filter for requested genes and sort by offset
            selected_guides = []
            processed_count = 0

            for gene in genes:
                if progress.wasCanceled():
                    return

                gene_guides = data[
                    (
                        data["locus_tag"].str.contains(gene)
                        | data["gene"].str.contains(gene)
                    )
                    & (data["mismatches"] == 0)
                ].sort_values("offset")

                # Take the specified number of guides
                selected_guides.append(gene_guides.head(self.guides_per_gene.value()))
                processed_count += self.guides_per_gene.value()
                progress.setValue(processed_count)

            if not selected_guides:
                raise ValueError("No matching guides found for the specified genes")

            # Combine all selected guides
            selected_guides = pd.concat(selected_guides)

            # Use temporary file for intermediate step
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=True
            ) as tmp_file:
                selected_guides.to_csv(tmp_file.name, sep="\t", index=False)

                # Configure and run mismatch.py with progress monitoring
                cmd = [
                    "python",
                    "mismatch.py",
                    "--parameters_file",
                    self.params_input.text(),
                    "--num_variants",
                    str(self.num_variants.value()),
                ]

                process = subprocess.Popen(
                    cmd,
                    stdin=open(tmp_file.name),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1,
                    universal_newlines=True,
                )

                # Buffer for collecting output
                output_buffer = StringIO()

                while True:
                    if progress.wasCanceled():
                        process.terminate()
                        return

                    line = process.stdout.readline()
                    if not line and process.poll() is not None:
                        break

                    output_buffer.write(line)
                    # Update progress based on processed variants
                    if line.strip():  # Only count non-empty lines
                        processed_count += 1
                        progress.setValue(min(processed_count, total_operations))

                # Check for errors
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(
                        process.returncode, cmd, process.stderr.read()
                    )

            # Store the results in memory
            output_buffer.seek(0)
            self.generated_mismatches = pd.read_csv(output_buffer, sep="\t")
            output_buffer.close()

            # Filter out the original guides
            self.generated_mismatches = self.generated_mismatches[
                self.generated_mismatches["mismatches"] > 0
            ]

            # Enable buttons and update UI with clearer titles
            self.save_btn.setEnabled(True)
            self.check_offtargets_btn.setEnabled(True)
            self.off_target_settings_box.setEnabled(True)
            self.genome_input.setEnabled(True)
            self.genome_input.setPlaceholderText(
                "Select genome file for off-target analysis"
            )
            self.browse_genome_btn.setEnabled(True)  # Updated reference
            self.off_target_settings_box.findChild(QLabel).setText(
                "Off-target Settings"
            )

            self.update_table_with_results()

            QMessageBox.information(
                self,
                "Success",
                f"Generated {len(self.generated_mismatches)} mismatched guides. Click Save to write to file.",
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
        finally:
            progress.close()

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

    def browse_genome_file(self):
        dialog = PreviewFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters(["GenBank Files (*.gb *.gbk)", "All Files (*)"])
        if dialog.exec_():
            chosen = dialog.selectedFiles()
            if chosen:
                self.genome_input.setText(chosen[0])
                self.check_offtargets_btn.setEnabled(True)

    def check_offtargets(self):
        """Check for off-targets in the generated mismatches"""
        if self.generated_mismatches is None:
            QMessageBox.warning(self, "No Data", "Generate mismatches first.")
            return

        genome_file = self.genome_input.text()
        if not genome_file:
            QMessageBox.warning(
                self, "No Genome File", "Please select a genome file first."
            )
            return

        try:
            # Create progress dialog before starting process
            self.progress = QProgressDialog(
                "Checking off-targets...",
                "Cancel",
                0,
                len(self.generated_mismatches),
                self,
            )
            self.progress.setWindowModality(Qt.WindowModal)
            self.progress.canceled.connect(self.cancel_process)
            self.progress.show()

            pam_direction = self.pam_direction_combo.currentText()

            # Create temporary FASTA file and store its path
            temp_fasta = tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            )
            self.temp_fasta_path = temp_fasta.name  # Store path before writing

            with temp_fasta:
                self.sequence_map = {}
                self.is_downstream = (
                    self.pam_direction_combo.currentText() == "downstream"
                )

                # print("\nDEBUG: Writing sequences to FASTA")
                for idx, row in self.generated_mismatches.iterrows():
                    spacer = row["spacer"]
                    pam = row["pam"]

                    if self.is_downstream:
                        # Store the original spacer keyed by the spacer sequence
                        self.sequence_map[spacer] = spacer
                        # For downstream PAMs, reverse complement the entire sequence for search
                        search_seq = self.rev_comp(spacer + pam)
                        # print(f"Original spacer={spacer} PAM={pam}")
                        # print(f"Search sequence: {search_seq}")
                    else:
                        # For upstream PAMs, use sequence as-is
                        search_seq = spacer
                        self.sequence_map[spacer] = spacer
                        # print(f"Original spacer={spacer} PAM={pam}")
                        # print(f"Search sequence: {search_seq}")

                    temp_fasta.write(f">{search_seq}\n{search_seq}\n")

            # Set up QProcess
            self.process = QProcess(self)
            self.process.setProcessChannelMode(QProcess.SeparateChannels)
            self.stdout_buffer = StringIO()
            self.process.readyReadStandardOutput.connect(self.handle_stdout)
            self.process.readyReadStandardError.connect(self.handle_stderr)
            self.process.finished.connect(self.offtarget_process_finished)

            # Run targets_spy.py with correct arguments
            mismatches = self.mismatches_combo.currentText()
            cmd = [
                sys.executable,
                "targets_spy.py",
                self.temp_fasta_path,  # Use stored path
                genome_file,
                "NGG",  # Always NGG since we handle orientation in sequence prep
                mismatches,
                "--pam-direction",
                "upstream",  # Always upstream since we handle orientation in sequence prep
            ]

            # print("\nDEBUG: Running command:", " ".join(cmd))
            self.process.start(cmd[0], cmd[1:])

        except Exception as e:
            self.cleanup()
            if hasattr(self, "progress"):
                self.progress.close()
            QMessageBox.critical(
                self, "Error", f"Failed to start off-target check:\n{str(e)}"
            )

    def handle_stdout(self):
        """Handle process stdout"""
        output = self.process.readAllStandardOutput().data().decode("utf-8")
        self.stdout_buffer.write(output)

    def handle_stderr(self):
        """Handle process stderr"""
        stderr = self.process.readAllStandardError().data().decode()
        if stderr:
            self.console.print(
                stderr, end=""
            )  # Using end="" to prevent double newlines

    def cancel_process(self):
        """Cancel the running process"""
        if self.process and self.process.state() != QProcess.NotRunning:
            self.process.kill()
        self.cleanup()

    def cleanup(self):
        """Clean up temporary files"""
        if hasattr(self, "temp_fasta_path") and os.path.exists(self.temp_fasta_path):
            try:
                os.unlink(self.temp_fasta_path)
            except:
                pass

    def offtarget_process_finished(self):
        """Handle completion of off-target checking"""
        self.progress.close()
        try:
            # Reset buffer position and read results
            self.stdout_buffer.seek(0)
            off_target_data = pd.read_csv(self.stdout_buffer, sep="\t")
            self.stdout_buffer.close()

            # print("\nDEBUG: Off-target analysis results")
            # print("=====================================")
            # print(f"Total guides analyzed: {len(off_target_data)}")
            # print("\nFirst few rows of off-target data:")
            # print(off_target_data.head())
            # print("\nGuide counts by number of sites:")
            # print(off_target_data["sites"].value_counts().sort_index())

            # Find guides with unique targets
            unique_targets = off_target_data[
                (off_target_data["sites"] == 1)
                & (off_target_data["mismatches"] == 1)
                & (off_target_data["intergenic"] == 0)
            ]

            # console log everything that is excluded
            print("Guides excluded:\n")
            excluded_guides = off_target_data[
                ~off_target_data["spacer"].isin(unique_targets["spacer"])
            ][
                [
                    "spacer",
                    "target",
                    "locus_tag",
                    "mismatches",
                    "sites",
                    "genes",
                    "intergenic",
                ]
            ].sort_values(
                "spacer"
            )  # Sort by spacer

            # Create rich table with box drawing characters for grouping
            table = Table(
                show_header=True,
                header_style="bold",
                box=DOUBLE_EDGE,  # Changed from box.DOUBLE_EDGE to DOUBLE_EDGE
                show_edge=True,  # Show outer border
                padding=(0, 1),  # Add some horizontal padding
            )

            for col in excluded_guides.columns:
                table.add_column(col)

            current_group = []
            current_spacer = None

            # Group rows by spacer
            for _, row in excluded_guides.iterrows():
                if current_spacer != row["spacer"]:
                    # Print previous group
                    if current_group:
                        for group_row in current_group:
                            table.add_row(*[str(x) for x in group_row])
                        table.add_section()  # Add a visual separator between groups
                    current_spacer = row["spacer"]
                    current_group = []
                current_group.append(row)

            # Add the last group
            if current_group:
                for group_row in current_group:
                    table.add_row(*[str(x) for x in group_row])

            console = Console()
            console.print(table)

            # print(f"\nGuides with unique targets: {len(unique_targets)}")
            # print("\nFirst few unique target sequences:")
            # print(unique_targets["spacer"].head())

            # Map the search sequences back to original spacers correctly
            matching_spacers = []
            # print("\nDEBUG: Processing unique targets:")
            for seq in unique_targets["spacer"]:
                if self.is_downstream:
                    # For downstream PAMs:
                    # 1. Rev comp the search result to get back to original orientation
                    orig_seq = self.rev_comp(seq)
                    # 2. Look up the original spacer
                    orig_spacer = self.sequence_map.get(orig_seq)
                    # print(f"Search result: {seq}")
                    # print(f"Rev comp back: {orig_seq}")
                    # print(f"Original spacer: {orig_spacer}")
                else:
                    # For upstream PAMs, use sequence as-is
                    orig_spacer = self.sequence_map.get(seq)
                    # print(f"Search result: {seq}")
                    # print(f"Original spacer: {orig_spacer}")

                if orig_spacer:
                    matching_spacers.append(orig_spacer)

            # Store counts for message
            initial_count = len(self.generated_mismatches)
            # print(f"\nInitial mismatches count: {initial_count}")

            # Debug guides being removed
            guides_to_remove = self.generated_mismatches[
                ~self.generated_mismatches["spacer"].isin(matching_spacers)
            ]
            # print("\nGuides being removed:")
            # print(
            #     guides_to_remove[
            #         [
            #             "spacer",
            #             "target",
            #             "pam",
            #             "locus_tag",
            #             "sites",
            #             "genes",
            #             "intergenic",
            #         ]
            #     ]
            # )

            # Filter using the original spacers
            self.generated_mismatches = self.generated_mismatches[
                self.generated_mismatches["spacer"].isin(matching_spacers)
            ]

            remaining_count = len(self.generated_mismatches)
            removed_count = initial_count - remaining_count

            # print(f"\nRemaining guides: {remaining_count}")
            # print(f"Removed guides: {removed_count}")

            # Update the table display
            self.update_table_with_results()

            QMessageBox.information(
                self,
                "Off-target Check Complete",
                f"Removed {removed_count} guides with off-targets.\n"
                f"Retained {remaining_count} guides with unique targets.\n\n"
                "Check the console for detailed debugging information.",
            )

        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to process off-target results:\n{str(e)}"
            )
        finally:
            self.cleanup()

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

    def load_genes_with_preview(self):
        dialog = PreviewFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters(["Text Files (*.txt)", "All Files (*)"])
        if dialog.exec_():
            chosen = dialog.selectedFiles()
            if chosen:
                file_path = chosen[0]
                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        file_content = f.read().strip()
                    self.genes_input.setPlainText(file_content)
                    self.genes_file_input.setText(
                        file_path
                    )  # Update the file path display
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"Failed to load file: {str(e)}")

    def rev_comp(self, sequence):
        """Return the reverse complement of a sequence"""
        return str(Seq(sequence).reverse_complement())


def main():
    app = QApplication(sys.argv)
    gui = MismatchDesignerGUI()
    gui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
