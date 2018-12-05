import logging
import os


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    def _comma_delimited_arg(self, key):
        if self._arguments[key]:
            return self._arguments[key].split(',')
        return []

    @property
    def in_dir(self):
        return os.path.expanduser(self._arguments['<in_dir>']) or None

    @property
    def out_dir(self):
        return os.path.expanduser(self._arguments['<out_dir>']) or None

    @property
    def exon_pos_file(self):
        return os.path.expanduser(self._arguments['<exon_pos_file>']) or None

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)
