import sys
import os
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QPushButton,
    QLabel,
    QStackedWidget,
    QStyle,
    QLineEdit,
    QMessageBox,
)
from PyQt5.QtGui import QIcon, QPixmap, QPalette, QColor
from PyQt5.QtCore import Qt
from targets_gui import BarcodeTargetSeekerGUI  # We'll move your existing GUI here
from find_guides_gui import FindGuidesGUI  # 1) Import
from assembly_finder_gui import AssemblyFinderGUI  # New import
from design_mismatches_gui import MismatchDesignerGUI  # New import

# Configure and start process
os.environ["COLUMNS"] = str(120)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Barcoder Suite")
        # self.setGeometry(100, 100, 600, 600)

        # Try to load icon if it exists
        icon_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "assets", "barcoder_icon.png"
        )
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))
        else:
            print(f"Icon not found at: {icon_path}")

        # Create stacked widget to handle different pages
        self.stacked_widget = QStackedWidget()
        self.setCentralWidget(self.stacked_widget)

        # Create and add pages
        self.welcome_page = self.create_welcome_page()
        self.targets_page = BarcodeTargetSeekerGUI()
        self.find_guides_page = FindGuidesGUI()  # 2) Instantiate
        self.assembly_page = AssemblyFinderGUI()  # Create new page
        self.mismatch_designer_page = MismatchDesignerGUI()  # Add new page instance

        self.stacked_widget.addWidget(self.welcome_page)
        self.stacked_widget.addWidget(self.targets_page)
        self.stacked_widget.addWidget(self.find_guides_page)
        self.stacked_widget.addWidget(self.assembly_page)
        self.stacked_widget.addWidget(self.mismatch_designer_page)  # Add to stack

        # Start with welcome page
        self.stacked_widget.setCurrentWidget(self.welcome_page)

    def create_welcome_page(self):
        page = QWidget()
        layout = QVBoxLayout()

        # Add Git Pull button at the top
        git_pull_btn = QPushButton("Update Software (Git Pull)")
        git_pull_btn.clicked.connect(self.perform_git_pull)
        git_pull_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #6c757d;
                color: white;
                padding: 10px;
                border-radius: 4px;
                font-weight: bold;
                max-width: 200px;
            }
            QPushButton:hover {
                background-color: #5a6268;
            }
        """
        )
        layout.addWidget(git_pull_btn, alignment=Qt.AlignRight)

        # Welcome message
        welcome_label = QLabel("Welcome to Barcoder Suite")
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; margin: 20px;")

        # Try to load and display logo
        logo_path = os.path.join(
            os.path.dirname(__file__), "assets", "barcoder_logo.png"
        )
        if os.path.exists(logo_path):
            logo_label = QLabel()
            pixmap = QPixmap(logo_path)
            scaled_pixmap = pixmap.scaled(
                200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation
            )
            logo_label.setPixmap(scaled_pixmap)
            logo_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(logo_label)

        # Load and display icon in the GUI
        icon_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "assets", "barcoder_icon.png"
        )
        if os.path.exists(icon_path):
            icon_label = QLabel()
            pixmap = QPixmap(icon_path)
            scaled_pixmap = pixmap.scaled(
                100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation
            )
            icon_label.setPixmap(scaled_pixmap)
            icon_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(icon_label)
        else:
            print(f"Icon not found at: {icon_path}")

        layout.addWidget(welcome_label)

        # Tool buttons
        find_guides_btn = self.create_tool_button(
            "Design New Guides",
            "Design new CRISPR-type barcodes.",
            lambda: self.stacked_widget.setCurrentWidget(self.find_guides_page),
        )
        layout.addWidget(find_guides_btn)

        targets_btn = self.create_tool_button(
            "Identify Guide Targets",
            "Find and analyze existing CRISPR-type barcodes.",
            lambda: self.stacked_widget.setCurrentWidget(self.targets_page),
        )

        layout.addWidget(targets_btn)

        mismatch_btn = self.create_tool_button(
            "Design Guide Mismatches",
            "Generate mismatched variants of existing CRISPR guides.",
            lambda: self.stacked_widget.setCurrentWidget(self.mismatch_designer_page),
        )
        layout.addWidget(mismatch_btn)

        assembly_btn = self.create_tool_button(
            "Assembly Finder",
            "Download sequences from NCBI by accession.",
            lambda: self.stacked_widget.setCurrentWidget(self.assembly_page),
        )
        layout.addWidget(assembly_btn)

        # Add stretch to push everything to the top
        layout.addStretch()

        page.setLayout(layout)
        return page

    def perform_git_pull(self):
        try:
            import subprocess

            result = subprocess.run(
                ["git", "pull"], capture_output=True, text=True, check=True
            )
            QMessageBox.information(
                self,
                "Git Pull Complete",
                f"Successfully updated software:\n{result.stdout}",
            )
            # Recommend restart
            restart = QMessageBox.question(
                self,
                "Restart Required",
                "Software has been updated. Would you like to restart the application?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if restart == QMessageBox.Yes:
                QApplication.quit()
                subprocess.Popen([sys.executable] + sys.argv)
        except subprocess.CalledProcessError as e:
            QMessageBox.critical(
                self, "Git Pull Failed", f"Failed to update software:\n{e.stderr}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred:\n{str(e)}")

    def create_tool_button(self, title, description, callback):
        button = QPushButton()
        button.clicked.connect(callback)

        # Create a container widget
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setSpacing(8)  # Increased spacing between elements
        layout.setContentsMargins(15, 15, 15, 15)  # Increased padding

        # Title with larger font and more spacing
        title_label = QLabel(title)
        title_label.setStyleSheet(
            """
            font-weight: bold;
            font-size: 18px;
            color: #212529;
            margin-bottom: 8px;
            background: transparent;
        """
        )

        # Description with better spacing
        desc_label = QLabel(description)
        desc_label.setStyleSheet(
            """
            color: #666;
            font-size: 13px;
            background: transparent;
            margin-top: 4px;
        """
        )
        desc_label.setWordWrap(True)

        layout.addWidget(title_label)
        layout.addWidget(desc_label)

        # Custom button with fixed size
        button = QPushButton()
        button.setFixedSize(400, 100)  # Set a comfortable size
        button.clicked.connect(callback)

        # Set the container as the button's content
        button.setLayout(layout)

        # Style the button
        button.setStyleSheet(
            """
            QPushButton {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 6px;
                text-align: left;
                padding: 15px;
            }
            QPushButton:hover {
                background-color: #e9ecef;
                border-color: #ced4da;
            }
        """
        )

        return button


def main():
    # Enable High DPI support
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor("#FFFFFF"))
    palette.setColor(QPalette.WindowText, QColor("#000000"))
    app.setPalette(palette)

    # Set application-wide stylesheet
    app.setStyleSheet(
        """
        QMainWindow { background-color: #ffffff; }
        QLineEdit, QComboBox, QSpinBox, QTextEdit, QPlainTextEdit {
            background-color: #ffffff;
            color: #000000;
            border: 1px solid #dee2e6;
        }
        QLabel { color: #212529; }
        QPushButton {
            background-color: #f8f9fa;
            color: #000000;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
        """
    )

    window = MainWindow()
    window.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
