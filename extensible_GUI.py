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
    QHBoxLayout,
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

        # Create top button row
        top_buttons = QHBoxLayout()

        # Relaunch button on the left
        relaunch_btn = QPushButton("Relaunch Application")
        relaunch_btn.clicked.connect(self.relaunch_app)
        relaunch_btn.setStyleSheet(
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
        top_buttons.addWidget(relaunch_btn, alignment=Qt.AlignLeft)

        # Add stretch to push buttons to opposite sides
        top_buttons.addStretch()

        # Git Pull button on the right
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
        top_buttons.addWidget(git_pull_btn, alignment=Qt.AlignRight)

        # Add the top button row to main layout
        layout.addLayout(top_buttons)

        # Create a container widget for buttons with center alignment
        button_container = QWidget()
        button_layout = QVBoxLayout(button_container)
        button_layout.setAlignment(Qt.AlignCenter)

        # Main Tools Section
        main_tools_label = QLabel("Main Tools")
        main_tools_label.setStyleSheet(
            "font-size: 18px; font-weight: bold; margin: 15px 0; color: #212529;"
        )
        main_tools_label.setAlignment(Qt.AlignCenter)
        button_layout.addWidget(main_tools_label)

        # Main tool buttons
        find_guides_btn = self.create_tool_button(
            "Design New Guides",
            "Design new CRISPR-type barcodes from scratch using your genome sequence.",
            lambda: self.stacked_widget.setCurrentWidget(self.find_guides_page),
        )
        button_layout.addWidget(find_guides_btn, 0, Qt.AlignCenter)

        mismatch_btn = self.create_tool_button(
            "Design Guide Mismatches",
            "Generate mismatched variants of existing CRISPR guides to create barcodes.",
            lambda: self.stacked_widget.setCurrentWidget(self.mismatch_designer_page),
        )
        button_layout.addWidget(mismatch_btn, 0, Qt.AlignCenter)

        # Extra Tools Section
        extra_tools_label = QLabel("Extra Tools")
        extra_tools_label.setStyleSheet(
            "font-size: 18px; font-weight: bold; margin: 15px 0; color: #666666;"
        )
        extra_tools_label.setAlignment(Qt.AlignCenter)
        button_layout.addWidget(extra_tools_label)

        targets_btn = self.create_tool_button(
            "Identify Guide Targets",
            "For use with existing guides: Analyze where guides from other sources might target in your genome.",
            lambda: self.stacked_widget.setCurrentWidget(self.targets_page),
        )
        button_layout.addWidget(targets_btn, 0, Qt.AlignCenter)

        assembly_btn = self.create_tool_button(
            "Assembly Finder",
            "Download sequences from NCBI by accession number.",
            lambda: self.stacked_widget.setCurrentWidget(self.assembly_page),
        )
        button_layout.addWidget(assembly_btn, 0, Qt.AlignCenter)

        # Add the button container to the main layout
        layout.addWidget(button_container)

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
        # Custom button with adjusted size
        button = QPushButton()
        button.setFixedWidth(400)  # Fixed width
        button.setMinimumHeight(120)  # Minimum height, will expand if needed
        button.clicked.connect(callback)

        # Create a container widget
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setSpacing(4)  # Reduced spacing
        layout.setContentsMargins(10, 10, 10, 10)  # Reduced margins

        # Title with larger font
        title_label = QLabel(title)
        title_label.setStyleSheet(
            """
            font-weight: bold;
            font-size: 16px;
            color: #212529;
            background: transparent;
        """
        )
        title_label.setWordWrap(True)

        # Description with proper wrapping
        desc_label = QLabel(description)
        desc_label.setStyleSheet(
            """
            color: #666;
            font-size: 12px;
            background: transparent;
        """
        )
        desc_label.setWordWrap(True)
        desc_label.setFixedWidth(360)  # Leave margin for padding

        layout.addWidget(title_label)
        layout.addWidget(desc_label)
        layout.addStretch()  # Push content to top

        # Set the container as the button's content
        button.setLayout(layout)

        # Update button style
        button.setStyleSheet(
            """
            QPushButton {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 6px;
                text-align: left;
                padding: 12px;
            }
            QPushButton:hover {
                background-color: #e9ecef;
                border-color: #ced4da;
            }
        """
        )

        return button

    def relaunch_app(self):
        """Relaunch the application"""
        try:
            import subprocess

            QApplication.quit()
            subprocess.Popen([sys.executable] + sys.argv)
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to relaunch application:\n{str(e)}"
            )
            # If relaunch fails, don't quit the current instance
            return


def main():
    # Force light mode based on platform
    if sys.platform.startswith("win"):
        sys.argv += ["-platform", "windows:darkmode=1"]
    elif sys.platform.startswith("darwin"):  # macOS
        sys.argv += ["-platform", "cocoa:darkmode=1"]
    else:  # Linux and others
        sys.argv += ["-platform", "xcb:darkmode=1"]

    # Enable High DPI support
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    # Force light theme palette
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor("#FFFFFF"))
    palette.setColor(QPalette.WindowText, QColor("#000000"))
    palette.setColor(QPalette.Base, QColor("#FFFFFF"))
    palette.setColor(QPalette.AlternateBase, QColor("#F5F5F5"))
    palette.setColor(QPalette.ToolTipBase, QColor("#FFFFFF"))
    palette.setColor(QPalette.ToolTipText, QColor("#000000"))
    palette.setColor(QPalette.Text, QColor("#000000"))
    palette.setColor(QPalette.Button, QColor("#F0F0F0"))
    palette.setColor(QPalette.ButtonText, QColor("#000000"))
    palette.setColor(QPalette.BrightText, QColor("#FF0000"))
    palette.setColor(QPalette.Link, QColor("#0000FF"))
    palette.setColor(QPalette.Highlight, QColor("#308CC6"))
    palette.setColor(QPalette.HighlightedText, QColor("#FFFFFF"))

    # Force the palette on the entire application
    app.setPalette(palette)

    # Set application-wide stylesheet to ensure light theme
    app.setStyleSheet(
        """
        QWidget {
            background-color: #ffffff;
            color: #000000;
        }
        QLineEdit, QComboBox, QSpinBox, QTextEdit, QPlainTextEdit {
            background-color: #ffffff;
            color: #000000;
            border: 1px solid #dee2e6;
        }
        QLabel { 
            color: #212529;
            background-color: transparent;
        }
        QPushButton {
            background-color: #f8f9fa;
            color: #000000;
        }
        QPushButton:hover {
            background-color: #e9ecef;
        }
        QTableWidget {
            background-color: #ffffff;
            color: #000000;
        }
        QHeaderView::section {
            background-color: #f8f9fa;
            color: #000000;
        }
        QComboBox QAbstractItemView {
            background-color: #ffffff;
            color: #000000;
        }
        QMessageBox {
            background-color: #ffffff;
        }
        QProgressDialog {
            background-color: #ffffff;
        }
    """
    )

    window = MainWindow()
    window.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
