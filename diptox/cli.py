# diptox/cli.py
import sys
import os
from pathlib import Path
from streamlit.web import cli as stcli


def _suppress_streamlit_welcome():
    try:
        streamlit_dir = Path.home() / ".streamlit"
        credentials_path = streamlit_dir / "credentials.toml"

        if not credentials_path.exists():
            streamlit_dir.mkdir(parents=True, exist_ok=True)
            with open(credentials_path, "w", encoding="utf-8") as f:
                f.write('[general]\nemail = ""\n')
    except Exception:
        pass


def run_gui():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    app_path = os.path.join(base_dir, "web_ui.py")

    if not os.path.exists(app_path):
        print(f"Error: Could not find GUI script at {app_path}")
        sys.exit(1)
    os.environ["DIPTOX_GUI_MODE"] = "true"
    _suppress_streamlit_welcome()
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