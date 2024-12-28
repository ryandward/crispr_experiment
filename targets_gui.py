import sys
import os
import argparse
import subprocess
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QComboBox,
    QMessageBox,
    QCheckBox,
    QProgressDialog,
    QFormLayout,
    QTextEdit,  # Added
)
from PyQt5.QtCore import Qt, QTimer, QRegExp, QProcess  # Added QProcess
from PyQt5.QtGui import QIntValidator, QRegExpValidator
import re
from rich.console import Console
from rich import print as rich_print


class BarcodeTargetSeekerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.parser = self.create_parser()
        self.process = None  # Added process attribute
        self.initUI()

    def create_parser(self):
        parser = argparse.ArgumentParser(
            description="Map barcodes to a circular genome"
        )
        parser.add_argument("sgrna_file", help="Path to sgrna_fasta_file", type=str)
        parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
        parser.add_argument("pam", help="PAM sequence", type=str)
        parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)
        parser.add_argument(
            "--pam-direction",
            choices=["upstream", "downstream"],
            default="downstream",
            help="Direction of the PAM sequence",
        )
        parser.add_argument(
            "--regulatory",
            type=int,
            default=0,
            help="Regulatory region size around target position.",
        )
        return parser

    def initUI(self):
        main_layout = QVBoxLayout()

        # Add a back button
        back_btn = QPushButton("← Back to Main Menu")
        back_btn.clicked.connect(
            lambda: self.parent().setCurrentWidget(
                self.parent().widget(0)  # Assuming welcome page is at index 0
            )
        )
        main_layout.addWidget(back_btn)

        # Title
        title = QLabel("Barcode Target Seeker")
        title.setStyleSheet("font-size: 20px; font-weight: bold; margin: 10px 0;")
        main_layout.addWidget(title)

        # Dynamic argument input fields
        self.argument_widgets = {}

        # All arguments
        all_actions = [
            action
            for action in self.parser._actions
            if action.dest not in ["help", None, "json"]
        ]

        # Use QFormLayout for inputs
        form_layout = QFormLayout()

        # Arguments
        for action in all_actions:
            if action.dest == "sgrna_file":
                label_text = "sgRNA File"
            else:
                label_text = action.dest.replace("_", " ").title()
            label = QLabel(label_text)
            # label.setMinimumWidth(120)

            if action.dest in ["sgrna_file", "genome_file"]:
                input_widget = QLineEdit()
                browse_button = QPushButton("Browse")
                browse_button.clicked.connect(
                    lambda checked, dest=action.dest: self.browse_file(dest)
                )
                file_layout = QHBoxLayout()
                file_layout.addWidget(input_widget)
                file_layout.addWidget(browse_button)
                form_layout.addRow(label, file_layout)
                self.argument_widgets[action.dest] = input_widget

            elif action.dest == "mismatches":
                input_widget = QComboBox()
                input_widget.addItems(["0", "1", "2"])
                input_widget.setCurrentText("1")  # Set default to 0
                form_layout.addRow(label, input_widget)
                self.argument_widgets[action.dest] = input_widget

            elif action.dest == "pam":
                input_widget = QLineEdit("NGG")  # Set default to NGG
                # Add validator for PAM sequence (only ACGTN letters)
                validator = QRegExpValidator(QRegExp("^[ACGTNacgtn]+$"))
                input_widget.setValidator(validator)
                form_layout.addRow(label, input_widget)
                self.argument_widgets[action.dest] = input_widget

            elif action.dest == "pam_direction":
                input_widget = QComboBox()
                input_widget.addItems(action.choices)
                input_widget.setCurrentText("downstream")  # default
                form_layout.addRow(label, input_widget)
                self.argument_widgets[action.dest] = input_widget

            elif action.dest == "regulatory":
                input_widget = QLineEdit("0")  # default as string
                form_layout.addRow(label, input_widget)
                self.argument_widgets[action.dest] = input_widget

        main_layout.addLayout(form_layout)

        # Add input field for output filename
        self.output_filename_label = QLabel("Output Filename:")
        self.output_filename_input = QLineEdit()
        self.output_filename_input.setPlaceholderText("Enter output filename")
        main_layout.addWidget(self.output_filename_label)
        main_layout.addWidget(self.output_filename_input)

        # Run button
        run_button = QPushButton("Run Barcode Target Seeker")
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

        # Add stretch to push everything to the top
        main_layout.addStretch()

        self.setLayout(main_layout)

    def browse_file(self, dest):
        if dest == "sgrna_file":
            filter_str = "FASTA Files (*.fasta *.fa);;All Files (*)"
        elif dest == "genome_file":
            filter_str = "GenBank Files (*.gb *.gbk);;All Files (*)"
        else:
            filter_str = "All Files (*)"
        filename, _ = QFileDialog.getOpenFileName(
            self,
            f"Select {self.argument_widgets[dest].placeholderText()}",
            "",
            filter_str,
        )
        if filename:
            self.argument_widgets[dest].setText(filename)

    def run_script(self):
        # First validate output filename as it's required
        output_filename = self.output_filename_input.text().strip()
        if not output_filename:
            QMessageBox.warning(self, "Input Error", "Output Filename is required.")
            return

        # Set appropriate extension based on output format
        if not output_filename.lower().endswith(".tsv"):
            output_filename += ".tsv"

        # Check if the output file already exists
        if os.path.exists(output_filename):
            reply = QMessageBox.question(
                self,
                "Overwrite Confirmation",
                f"The file '{output_filename}' already exists. Do you want to overwrite it?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply == QMessageBox.No:
                QMessageBox.information(
                    self,
                    "Operation Cancelled",
                    "Please choose a different output filename.",
                )
                return
            else:
                try:
                    # Overwrite the existing file by opening it in write mode and closing immediately
                    with open(output_filename, "w") as f:
                        pass
                except Exception as e:
                    QMessageBox.critical(
                        self,
                        "Error",
                        f"Failed to overwrite the file:\n{str(e)}",
                    )
                    return

        # Validate PAM sequence
        pam_input = self.argument_widgets["pam"].text().upper().strip()
        if not pam_input or not re.match("^[ACGTN]+$", pam_input):
            QMessageBox.warning(
                self, "Input Error", "Invalid PAM sequence. Use A, C, G, T, or N."
            )
            return

        # Collect required positional arguments in the correct order
        sgrna_file = self.argument_widgets.get("sgrna_file").text().strip()
        genome_file = self.argument_widgets.get("genome_file").text().strip()
        mismatches = self.argument_widgets.get("mismatches").currentText().strip()

        if not sgrna_file or not genome_file or mismatches == "":
            QMessageBox.warning(
                self,
                "Input Error",
                "sgRNA File, Genome File, and Mismatches are required.",
            )
            return

        args = [
            sgrna_file,
            genome_file,
            pam_input,
            mismatches,
        ]

        pam_dir = self.argument_widgets["pam_direction"].currentText().strip()

        # Confirm direction if PAM ends with 'N' but user selected 'downstream'
        if pam_input.endswith("N") and pam_dir == "downstream":
            msg = (
                "Does your PAM really end with 'N' and is downstream from the protospacer?\n"
                "Click OK to continue or Cancel to revise."
            )
            reply = QMessageBox.question(
                self, "Confirm PAM Direction", msg, QMessageBox.Ok | QMessageBox.Cancel
            )
            if reply == QMessageBox.Cancel:
                return

        # Confirm direction if PAM starts with 'N' but user selected 'upstream'
        if pam_input.startswith("N") and pam_dir == "upstream":
            msg = (
                "Does your PAM really start with 'N' and is upstream from the protospacer?\n"
                "Click OK to continue or Cancel to revise."
            )
            reply = QMessageBox.question(
                self, "Confirm PAM Direction", msg, QMessageBox.Ok | QMessageBox.Cancel
            )
            if reply == QMessageBox.Cancel:
                return

        if pam_dir != "downstream":
            args.extend(["--pam-direction", pam_dir])

        reg_val = self.argument_widgets["regulatory"].text().strip()
        if reg_val != "0":
            args.extend(["--regulatory", reg_val])

        # Set up and start progress dialog
        self.progress = QProgressDialog(
            "Running Barcode Target Seeker...", "Cancel", 0, 0, self
        )
        self.progress.setWindowModality(Qt.WindowModal)
        self.progress.setAutoReset(False)
        self.progress.setAutoClose(False)
        self.progress.canceled.connect(self.cancel_process)
        self.progress.show()

        self.process = QProcess()
        self.process.setProcessChannelMode(QProcess.SeparateChannels)
        self.process.readyReadStandardOutput.connect(
            lambda: self.handle_stdout(output_filename)
        )
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(
            lambda exit_code, exit_status: self.process_finished(
                exit_code, output_filename
            )
        )

        script_path = os.path.join(os.path.dirname(__file__), "targets.py")
        self.process.start(sys.executable, [script_path] + args)

    def handle_stdout(self, output_filename):
        stdout = self.process.readAllStandardOutput().data().decode()
        if stdout:
            # Open the file in append mode to add new content without overwriting
            with open(output_filename, "a") as f:
                f.write(stdout)

    def handle_stderr(self):
        stderr = self.process.readAllStandardError().data().decode()
        if stderr:
            # Create a new console for each output to ensure proper formatting
            temp_console = Console(
                force_terminal=True,
                color_system="truecolor",
                width=None,
                soft_wrap=True,
                highlight=True,
            )
            temp_console.print(stderr.rstrip())

    def cancel_process(self):
        if self.process and self.process.state() != QProcess.NotRunning:
            self.process.kill()

    def process_finished(self, exit_code, output_filename):
        self.progress.close()

        if exit_code == 0:
            QMessageBox.information(
                self,
                "Success",
                f"Script executed successfully!\nOutput saved to: {output_filename}",
            )
        else:
            QMessageBox.critical(
                self, "Error", f"Script failed with exit code {exit_code}"
            )
