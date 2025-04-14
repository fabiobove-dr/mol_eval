import logging


class Logger:
    _instance = None

    def __new__(cls, log_level: int = logging.INFO):
        if cls._instance is None:
            cls._instance = super(Logger, cls).__new__(cls)
            cls._instance._initialize(log_level)
        return cls._instance

    def _initialize(self, log_level: int):
        """Initialize the logger."""
        self.logger = logging.getLogger("ProjectLogger")
        self.logger.setLevel(log_level)

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))

        # Add the console handler to the logger
        self.logger.addHandler(console_handler)

    def get_logger(self):
        """Return the logger instance."""
        return self.logger
