# diptox/user_reg.py
import json
import requests
import logging
from pathlib import Path


# Google Forms
GOOGLE_CONFIG = {
    "url": "https://docs.google.com/forms/d/e/1FAIpQLSeGcDAqyM0SV2bPwcRupBa5Jg_SLHWSmjYdERkd1EjxqqnkDA/formResponse",
    "fields": {
        "name": "entry.1009922341",
        "affiliation": "entry.183441193",
        "email": "entry.1851128859"
    }
}

CN_CONFIG = {
    "type": "feishu_webhook",
    "url": "https://open.feishu.cn/open-apis/bot/v2/hook/84d3dbd8-3264-4dd6-917a-e0a5bd111787",
}

CONFIG_DIR = Path.home() / ".diptox"
CONFIG_FILE = CONFIG_DIR / "user_info.json"


def is_google_accessible(timeout=3):
    try:
        requests.head("https://www.google.com", timeout=timeout)
        return True
    except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
        return False
    except Exception:
        return False


def is_registered_or_skipped():
    if not CONFIG_FILE.exists():
        return False
    try:
        with open(CONFIG_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
            return data.get("status") in ["registered", "skipped"]
    except:
        return False


def save_status(status, info=None):
    try:
        CONFIG_DIR.mkdir(parents=True, exist_ok=True)
        data = {"status": status, "info": info or {}}
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)
    except Exception as e:
        logging.error(f"Unable to save user status: {e}")


def submit_to_google(name, affiliation, email):
    data = {
        GOOGLE_CONFIG["fields"]["name"]: name,
        GOOGLE_CONFIG["fields"]["affiliation"]: affiliation,
        GOOGLE_CONFIG["fields"]["email"]: email
    }
    try:
        resp = requests.post(GOOGLE_CONFIG["url"], data=data, timeout=10)
        return resp.status_code == 200
    except:
        return False


def submit_to_feishu(name, affiliation, email):
    payload = {
        "msg_type": "text",
        "content": {
            "text": f"[DiPTox User Registration]\nName: {name}\nAffiliation: {affiliation}\nEmail: {email}"
        }
    }
    try:
        resp = requests.post(CN_CONFIG["url"], json=payload, timeout=10)
        return resp.status_code == 200
    except Exception as e:
        return False


def submit_info(name, affiliation, email):
    use_google = is_google_accessible()

    success = False
    msg = ""

    if use_google:
        logging.info("Detected Google access, utilizing Google Forms.")
        if submit_to_google(name, affiliation, email):
            success = True
            msg = "Submitted successfully via Google Forms."
        else:
            if submit_to_feishu(name, affiliation, email):
                success = True
                msg = "Google unavailable, submitted via backup channel."
    else:
        logging.info("Google unreachable, utilizing domestic channel.")
        if submit_to_feishu(name, affiliation, email):
            success = True
            msg = "Submitted successfully."

    if success:
        save_status("registered", {"name": name, "affiliation": affiliation, "email": email})
        return True, "Thank you! Information received."
    else:
        return False, "Connection failed. Please check your internet."
