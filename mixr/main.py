"""
MIXR: Mismatching Isoform eXon Remover

Usage:
  mixr <in_dir> <out_dir> <exon_pos_file> [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
import logging
import os
from mixr.constants import VERSION
from mixr.config import CommandLineArguments
from mixr.filterexons import filter_exons
from docopt import docopt


def main(**kwargs):
    docopt_args = docopt(__doc__, version=VERSION)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    filter_exons(arguments)


if __name__ == '__main__':
    main()
