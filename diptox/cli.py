# diptox/cli.py
import sys
import os
from streamlit.web import cli as stcli


def run_gui():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    app_path = os.path.join(base_dir, "web_ui.py")

    if not os.path.exists(app_path):
        print(f"Error: Could not find GUI script at {app_path}")
        sys.exit(1)
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--browser.gatherUsageStats", "false",
        "--server.headless", "false"
    ]
    sys.exit(stcli.main())


if __name__ == "__main__":
    run_gui()