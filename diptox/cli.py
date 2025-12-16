# diptox/cli.py
import sys
import os
import time
import platform
import errno
from pathlib import Path
from streamlit.web import cli as stcli

CONFIG_DIR = Path.home() / ".diptox"
TEMP_DIR = CONFIG_DIR / "temp"


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


def _is_pid_running(pid: int) -> bool:
    if platform.system() == "Windows":
        try:
            os.kill(pid, 0)
            return True
        except OSError as e:
            if getattr(e, 'winerror', None) == 87:
                return False
            if getattr(e, 'winerror', None) == 5:
                return True
            return True
        except Exception:
            return True
    else:
        # Unix/Linux/Mac
        try:
            os.kill(pid, 0)
            return True
        except OSError as e:
            if e.errno == errno.ESRCH:  # No such process
                return False
            if e.errno == errno.EPERM:  # Permission denied
                return True
            return True
        except Exception:
            return True


def _clean_residual_temp_files():
    if not TEMP_DIR.exists():
        return

    print(f"Checking for residual temp files in {TEMP_DIR}...")
    cleaned_count = 0
    current_pid = os.getpid()

    for item in TEMP_DIR.iterdir():
        try:
            if not item.is_file():
                continue

            name = item.name
            parts = name.split('_')

            should_delete = False
            if len(parts) >= 2 and parts[0] == 'tmp' and parts[1].isdigit():
                file_pid = int(parts[1])
                if file_pid == current_pid:
                    continue
                if _is_pid_running(file_pid):
                    continue
                should_delete = True
            elif time.time() - item.stat().st_mtime > 86400:
                should_delete = True

            if should_delete:
                try:
                    item.unlink()
                    cleaned_count += 1
                except PermissionError:
                    print(f"Skipping {name}: File is locked by another process.")
                except Exception as e:
                    print(f"Failed to delete {name}: {e}")

        except Exception:
            pass

    if cleaned_count > 0:
        print(f"Cleaned {cleaned_count} residual temporary file(s).")


def run_gui():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    app_path = os.path.join(base_dir, "web_ui.py")

    if not os.path.exists(app_path):
        print(f"Error: Could not find GUI script at {app_path}")
        sys.exit(1)
    os.environ["DIPTOX_GUI_MODE"] = "true"
    _suppress_streamlit_welcome()
    _clean_residual_temp_files()
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