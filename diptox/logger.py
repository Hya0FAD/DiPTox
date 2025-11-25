# diptox/logger.py
import logging
import sys
import time
from logging.handlers import RotatingFileHandler, TimedRotatingFileHandler
from pathlib import Path
from typing import Optional, Union
import re
import os

# Default log configuration
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
DEFAULT_LOG_DIR = "DiPTox-Logs"
DEFAULT_LOG_LEVEL = logging.INFO


class LogManager:
    def __init__(self):
        self._configured = False
        self.loggers = {}

    def configure(
        self,
        log_dir: Optional[Union[str, Path]] = None,
        console_level: int = logging.INFO,
        file_level: int = logging.DEBUG,
        max_bytes: int = 10 * 1024 * 1024,  # 10MB
        backup_count: int = 5,
        when: str = "midnight",
        interval: int = 1,
        fmt: str = DEFAULT_LOG_FORMAT,
        datefmt: str = DEFAULT_DATE_FORMAT,
        enable_console: bool = True,
        enable_file: bool = True,
        log_retention_days: int = 2,
        max_total_logs: Optional[int] = 5
    ):
        """Global log configuration method."""
        root_logger = logging.getLogger()

        if root_logger.hasHandlers():
            root_logger.handlers.clear()

        if os.environ.get("DIPTOX_GUI_MODE") == "true":
            enable_file = False

        log_dir = Path(log_dir) if log_dir else Path(DEFAULT_LOG_DIR)
        if enable_file:
            log_dir.mkdir(parents=True, exist_ok=True)

        # Basic log formatter
        formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)

        # Initialize the root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)

        # Console handler
        if enable_console:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(console_level)
            console_handler.setFormatter(formatter)
            root_logger.addHandler(console_handler)

        # File handler (rotating by size)
        if enable_file:
            file_handler = RotatingFileHandler(
                filename=log_dir / "DiPTox.log",
                maxBytes=max_bytes,
                backupCount=backup_count,
                encoding="utf-8"
            )
            file_handler.setLevel(file_level)
            file_handler.setFormatter(formatter)
            root_logger.addHandler(file_handler)

            # Separate handler for error logs (rotating by time)
            error_handler = TimedRotatingFileHandler(
                filename=log_dir / "DiPToxError.log",
                when=when,
                interval=interval,
                backupCount=backup_count,
                encoding="utf-8"
            )
            error_handler.setLevel(logging.WARNING)
            error_handler.setFormatter(formatter)
            root_logger.addHandler(error_handler)

            self._setup_log_cleaner(
                log_dir=log_dir,
                retention_days=log_retention_days,
                max_total=max_total_logs
            )

        # Configure log level for third-party libraries
        logging.getLogger("urllib3").setLevel(logging.WARNING)
        logging.getLogger("requests").setLevel(logging.WARNING)

        self._configured = True

    def get_logger(self, name: Optional[str] = None) -> logging.Logger:
        """Retrieve a configured logger."""
        if not self._configured:
            self.configure()

        if not name:
            name = "root"

        if name in self.loggers:
            return self.loggers[name]

        logger = logging.getLogger(name)
        self.loggers[name] = logger
        return logger

    def _setup_log_cleaner(self, log_dir: Path,
                           retention_days: int,
                           max_total: Optional[int]):
        """Configuring log clearing with backup pattern matching"""
        now = time.time()

        main_logs = {'DiPTox.log', 'DiPToxError.log'}
        backup_pattern = re.compile(r'^(.*?)\.(log\.\d+|log\.\d{4}-\d{2}-\d{2})$')

        backup_files = []
        for f in log_dir.glob("*"):
            if f.name in main_logs:
                continue
            if backup_pattern.match(f.name):
                backup_files.append(f)

        for f in backup_files:
            if f.stat().st_mtime < now - retention_days * 86400:
                try:
                    f.unlink()
                    logging.debug(f"Deleted old log backup: {f}")
                except PermissionError as e:
                    logging.warning(f"The log backup cannot be deleted {f}: {e}")

        if max_total and len(backup_files) > max_total:
            sorted_files = sorted(backup_files, key=lambda x: x.stat().st_mtime)
            for f in sorted_files[:len(backup_files) - max_total]:
                try:
                    f.unlink()
                    logging.debug(f"Removed excess log backup: {f}")
                except PermissionError as e:
                    logging.warning(f"The log backup cannot be deleted {f}: {e}")

        # Configure periodically clearing threads
        if not hasattr(self, '_cleaner_thread'):
            self._start_cleaner_thread(log_dir, retention_days)

    def _start_cleaner_thread(self, log_dir: Path, interval_days: int = 1):
        """Start the background cleanup thread"""
        from threading import Thread

        def cleaner():
            while True:
                time.sleep(interval_days * 86400)
                self._setup_log_cleaner(log_dir, interval_days, None)
        self._cleaner_thread = Thread(
            target=cleaner,
            name="LogCleaner",
            daemon=True
        )
        self._cleaner_thread.start()


# Singleton instance
log_manager = LogManager()
