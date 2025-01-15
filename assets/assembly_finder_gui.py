# assembly_finder_gui.py

import os
import sys
import re
import configparser
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QFormLayout,
    QLabel,
    QLineEdit,
    QComboBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QProgressDialog,
    QMessageBox,
    QTextEdit,
)
from PyQt5.QtCore import Qt, QProcess, QRegExp
from PyQt5.QtGui import QColor, QFont, QTextCursor


class AssemblyFinderGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.process = None
        self.stdout_buffer = ""
        self.stderr_buffer = ""
        self.current_data = []  # Store current data for table
        self.log_file = "download_log.tsv"  # Add log file path
        self.initUI()
        self.read_config()
        self.load_existing_log()  # Add loading of existing log

    def initUI(self):
        self.setWindowTitle("Assembly Finder")
        self.setGeometry(100, 100, 900, 700)

        main_layout = QVBoxLayout()

        # Back button at the top
        back_btn = QPushButton("← Back to Main Menu")
        back_btn.clicked.connect(
            lambda: self.parent().setCurrentWidget(
                self.parent().widget(0)  # Assuming welcome page is at index 0
            )
        )
        back_btn.setStyleSheet(
            """
            QPushButton {
                padding: 5px;
                margin-bottom: 10px;
            }
            """
        )
        main_layout.addWidget(back_btn)

        # Title
        title_label = QLabel("Assembly Finder")
        title_label.setStyleSheet("font-size: 20px; font-weight: bold; margin: 10px 0;")
        main_layout.addWidget(title_label)

        # Form layout for inputs
        form_layout = QFormLayout()

        self.accessions_edit = QLineEdit()
        form_layout.addRow("Accession Numbers (space separated):", self.accessions_edit)

        self.source_combo = QComboBox()
        self.source_combo.addItems(["GenBank", "RefSeq"])
        form_layout.addRow("Source:", self.source_combo)

        self.email_edit = QLineEdit()
        form_layout.addRow("Entrez Email:", self.email_edit)

        main_layout.addLayout(form_layout)

        # Run button with matching style
        run_button = QPushButton("Download Assemblies")
        run_button.setStyleSheet(
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
        run_button.clicked.connect(self.run_script)
        main_layout.addWidget(run_button)

        # Add cleanup button after the Run button
        cleanup_button = QPushButton("Clean Up Log")
        cleanup_button.setStyleSheet(
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
            """
        )
        cleanup_button.clicked.connect(self.cleanup_log)
        main_layout.addWidget(cleanup_button)

        # Table for results
        table_label = QLabel("Downloaded Assemblies:")
        main_layout.addWidget(table_label)

        self.table = QTableWidget()
        self.table.setColumnCount(4)  # Changed from 3 to 4
        self.table.setHorizontalHeaderLabels(
            ["Accession", "Organism", "File Path", "Source"]
        )
        self.table.horizontalHeader().setStretchLastSection(True)
        main_layout.addWidget(self.table)

        self.setLayout(main_layout)

    def read_config(self):
        self.config_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "config.ini"
        )
        self.config = configparser.ConfigParser()
        if os.path.exists(self.config_path):
            self.config.read(self.config_path)
        if not self.config.has_section("Entrez"):
            self.config.add_section("Entrez")
        saved_email = self.config.get("Entrez", "email", fallback="")
        self.email_edit.setText(saved_email)

    def write_config(self, email):
        self.config.set("Entrez", "email", email)
        with open(self.config_path, "w") as f:
            self.config.write(f)

    def run_script(self):
        accessions = self.accessions_edit.text().strip()
        if not accessions:
            QMessageBox.warning(
                self, "Input Error", "Please provide at least one accession number."
            )
            return

        email_arg = self.email_edit.text().strip()
        if not email_arg:
            QMessageBox.warning(
                self, "Input Error", "Please provide an Entrez email address."
            )
            return

        # Save email to config
        self.write_config(email_arg)

        source = self.source_combo.currentText()

        # Prepare arguments
        script_path = os.path.join(os.path.dirname(__file__), "assembly_finder.py")
        args = [
            "--gui",  # Enable simple output mode
            *accessions.split(),
            "--source",
            source,
            "--email",
            email_arg,
        ]

        # Initialize progress dialog
        self.progress = QProgressDialog("Initializing...", "Cancel", 0, 0, self)
        self.progress.setWindowModality(Qt.WindowModal)
        self.progress.setMinimumDuration(0)
        self.progress.canceled.connect(self.cancel_process)
        self.progress.show()

        self.progress.setLabelText("Starting download...")

        # Initialize QProcess
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.SeparateChannels)
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.process_finished)

        # Start the process
        self.process.start(sys.executable, [script_path] + args)

        if not self.process.waitForStarted(1000):
            QMessageBox.critical(
                self, "Error", "Failed to start the assembly_finder script."
            )
            self.progress.close()

    def handle_stdout(self):
        output = self.process.readAllStandardOutput().data().decode("utf-8")
        self.stdout_buffer += output
        # Currently, assembly_finder.py directs rich to stderr, so stdout may not be used.
        # If needed, handle stdout here.

    def handle_stderr(self):
        error_output = self.process.readAllStandardError().data().decode("utf-8")
        for line in error_output.splitlines():
            # Handle progress updates
            if "Fetching" in line:
                self.progress.setLabelText("Fetching assembly accessions...")
            elif "Downloading" in line:
                self.progress.setLabelText("Downloading sequences...")
            elif "Saving" in line:
                self.progress.setLabelText("Saving sequences...")
            elif "Download complete" in line:
                self.progress.setLabelText("Download complete!")

            # Store structured data for table
            if "\t" in line and not line.startswith("["):
                self.stderr_buffer += line + "\n"

    def strip_ansi_codes(self, text):
        ansi_escape = re.compile(r"(?:\x1B[@-_][0-?]*[ -/]*[@-~])")
        return ansi_escape.sub("", text)

    def update_progress_dialog(self, message):
        """
        Update the progress dialog's label based on the message content.
        """
        # Define patterns to look for in the stderr messages
        patterns = {
            "Fetching": r"\[1/3\]\s+Fetching assembly accessions...",
            "Downloading": r"\[2/3\]\s+Downloading sequences...",
            "Saving": r"\[3/3\]\s+Saving sequences to files...",
            "Complete": r"Download complete!",
        }

        # Check which pattern matches the current message
        for key, pattern in patterns.items():
            if re.search(pattern, message):
                if key == "Complete":
                    self.progress.setLabelText("Download Complete!")
                    self.progress.setValue(100)
                else:
                    # Extract the task description without the progress bar
                    # Using regex groups for better extraction
                    match = re.search(r"\[.\]/3\]\s+(.*)", message)
                    if match:
                        task_description = match.group(1).replace("...", "")
                        self.progress.setLabelText(task_description)
                break

    def load_existing_log(self):
        """Load existing download log into the table"""
        if not os.path.exists(self.log_file):
            return

        with open(self.log_file, "r") as f:
            lines = f.readlines()
            if not lines:
                return

            # Get header and data
            header = lines[0].strip().split("\t")[1:]  # Skip timestamp column
            data = []

            for line in lines[1:]:  # Skip header
                parts = line.strip().split("\t")[1:]  # Skip timestamp
                if len(parts) == 3:  # Old format without source
                    accession, organism, filepath = parts
                    # Determine source from accession prefix
                    source = (
                        "GenBank"
                        if accession.startswith("GCA")
                        else "RefSeq" if accession.startswith("GCF") else "Unknown"
                    )
                else:  # New format with source
                    accession, organism, filepath, source = parts

                # Convert absolute path to relative display path for table only
                display_path = os.path.join("./sequences", os.path.basename(filepath))
                data.append([accession, organism, display_path, source])

            # Set up table with 4 columns
            self.table.clear()
            self.table.setRowCount(len(data))
            self.table.setColumnCount(4)
            self.table.setHorizontalHeaderLabels(
                ["Accession", "Organism", "File Path", "Source"]
            )

            # Add data with alternating colors
            for row, row_data in enumerate(data):
                bg_color = QColor("#f8f9fa") if row % 2 == 0 else QColor("#ffffff")
                for col, value in enumerate(row_data):
                    item = QTableWidgetItem(value)
                    item.setBackground(bg_color)
                    self.table.setItem(row, col, item)

            # Adjust table appearance
            self.table.resizeColumnsToContents()
            self.table.horizontalHeader().setStretchLastSection(True)

    def cleanup_log(self):
        """Clean up log by removing duplicates and non-existent files"""
        if not os.path.exists(self.log_file):
            return

        # Read existing log
        with open(self.log_file, "r") as f:
            lines = f.readlines()
            if not lines:
                return

        # Process lines and keep track of unique entries
        header = lines[0]
        entries = {}  # Using dict to maintain latest entry for each accession

        for line in lines[1:]:  # Skip header
            parts = line.strip().split("\t")
            if len(parts) == 4:  # Old format without source
                timestamp, accession, organism, filepath = parts
                # Determine source from accession prefix
                source = (
                    "GenBank"
                    if accession.startswith("GCA")
                    else "RefSeq" if accession.startswith("GCF") else "Unknown"
                )
            else:  # New format with source
                timestamp, accession, organism, filepath, source = parts

            if os.path.exists(filepath):  # Only keep if file exists
                entries[accession] = (timestamp, organism, filepath, source)

        # Write cleaned log
        with open(self.log_file, "w") as f:
            f.write("Timestamp\tAccession\tOrganism\tFilePath\tSource\n")
            for accession, (timestamp, organism, filepath, source) in entries.items():
                f.write(f"{timestamp}\t{accession}\t{organism}\t{filepath}\t{source}\n")

        # Reload the table
        self.load_existing_log()

        # Show summary
        removed = len(lines) - len(entries) - 1  # -1 for header
        QMessageBox.information(
            self,
            "Log Cleaned",
            f"Removed {removed} entries from log:\n"
            f"- Deleted duplicates\n"
            f"- Removed entries for non-existent files\n"
            f"Current log entries: {len(entries)}",
        )

    def process_finished(self):
        self.progress.close()

        # Parse the tab-separated data from current run
        lines = self.stderr_buffer.strip().splitlines()
        header = None
        data = []

        for line in lines:
            if "AssemblyAccession" in line:
                header = line.split("\t")
            elif "\t" in line and not "Download complete" in line:
                data.append(line.split("\t"))

        if header and data:
            # Reload entire log including new data
            self.load_existing_log()

        QMessageBox.information(
            self, "Success", f"Successfully downloaded {len(data)} sequences."
        )

    def cancel_process(self):
        if self.process and self.process.state() != QProcess.NotRunning:
            self.process.kill()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    gui = AssemblyFinderGUI()
    gui.show()
    sys.exit(app.exec_())
