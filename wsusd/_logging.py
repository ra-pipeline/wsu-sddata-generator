import logging

LOG_FORMAT = '%(levelname)s: %(message)s'


def configure_logging(level=logging.INFO):
    if not logging.getLogger().handlers:
        logging.basicConfig(format=LOG_FORMAT, level=level)
    else:
        logging.getLogger().setLevel(level)


def get_logger(name):
    return logging.getLogger(name)


def set_debug_level(logger):
    logger.setLevel(logging.DEBUG)

