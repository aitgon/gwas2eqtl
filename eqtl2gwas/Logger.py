import logging

from termcolor import colored


class Logger:
    """
    This class defines the logger.

    """

    log_level = logging.DEBUG  # default log level

    logger = logging.getLogger("eqtl2gwas")
    formatter_str = '%(asctime)s :: %(levelname)s :: %(name)s :: %(message)s'
    handler = logging.StreamHandler()
    __formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: %(message)s')
    handler.setFormatter(__formatter)
    logger.addHandler(handler)
    logger.setLevel(log_level)  # set root's level

    @classmethod
    def debug(cls, msg):

        formatter_stream = logging.Formatter(colored(cls.formatter_str, 'blue', attrs=['bold']))
        cls.handler.setFormatter(formatter_stream)
        cls.logger.debug(msg)

    @classmethod
    def info(cls, msg):

        formatter_stream = logging.Formatter(colored(cls.formatter_str, 'green', attrs=['dark']))
        cls.handler.setFormatter(formatter_stream)
        cls.logger.info(msg)

    @classmethod
    def warning(cls, msg):

        formatter_stream = logging.Formatter(colored(cls.formatter_str, 'yellow', attrs=['bold']))
        cls.handler.setFormatter(formatter_stream)
        cls.logger.warning(msg)

    @classmethod
    def error(cls, msg):

        formatter_stream = logging.Formatter(colored(cls.formatter_str, 'red', attrs=[]))
        cls.handler.setFormatter(formatter_stream)
        cls.logger.error(msg)

    @classmethod
    def critical(cls, msg):

        formatter_stream = logging.Formatter(colored(cls.formatter_str, 'red', attrs=[]))
        cls.handler.setFormatter(formatter_stream)
        cls.logger.critical(msg)
