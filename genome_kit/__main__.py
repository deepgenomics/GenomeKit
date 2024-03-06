# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
#
# This file allows the genome_kit package to be run as a script.
# The script provides subcommands related to GenomeKit development
# and maintenance, such as building internal data files.
#
from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

from . import _build_annotations, _build_appris, gk_data

#################################################################


def parse_args():
    """Parses the command-line args when genome_kit is run as a script
    via ``python -m genome_kit ...``
    """
    # Create top-level parser and subcommand parsers
    parser = argparse.ArgumentParser(
        prog="genome_kit",
        description="""
        Provides subcommands related to GenomeKit development
        and maintenance, such as building internal data files.
        """.strip())

    subparsers = parser.add_subparsers(
        title="subcommands", help="use <command> -h for subcommand help", dest="subcommand")

    # Parser for "build" command
    build_parser = subparsers.add_parser(
        "build",
        description="The build subcommand builds full- or test-sized"
        "data files. Full-sized files can optionally be"
        "uploaded to the GenomeKit store. Test-sized"
        "files will be placed in the GenomeKit source tree"
        "under tests/data.",
        help="build data files, optionally uploading to store")
    build_parser.add_argument("--all", action="store_true", default=False, help="build all data files")
    build_parser.add_argument("--appris", action="store_true", default=False, help="build full-sized APPRIS data files")
    build_parser.add_argument(
        "--test-anno", action="store_true", default=False, help="build test-sized annotation data files")
    build_parser.add_argument(
        "--test-appris", action="store_true", default=False, help="build test-sized APPRIS data files")
    build_parser.add_argument(
        "--test-2bit", action="store_true", default=False, help="build test-sized 2bit data files")
    build_parser.add_argument(
        "--upload",
        action="store_true",
        default=False,
        help="upload any full-sized generated to GenomeKit store")

    # Parse all arguments and return the resulting namespace
    return parser.parse_args()


def build(args):
    """Builds GenomeKit's standard internal data files.

    Users can specify which subset of files they want built, and whether they want
    to upload the resulting files to the store.
    """

    # Build full-sized APPRIS
    if args.appris or args.all:
        _build_appris.build_full_appris_files(args.upload)

    _get_file_original = gk_data.data_manager.get_file
    test_data_dir = Path(__file__).parent.parent / 'tests' / 'data' / 'mini1'

    def get_test_file(filename):
        if (test_data_dir / filename).exists():
            return str(test_data_dir / filename)
        return _get_file_original(filename)

    gk_data.data_manager.get_file = get_test_file

    # Build test-sized 2bit files
    if args.test_2bit or args.all:
        _build_annotations.build_test_2bit_files()

    # Build test-sized annotation files
    if args.test_anno or args.all:
        _build_annotations.build_test_annotation_files()

    # Build test-sized annotation files
    if args.test_appris or args.all:
        _build_appris.build_test_appris_files()


def main():
    # Suppress progress messages from C++, e.g. when building annotations.
    os.environ["GENOMEKIT_QUIET"] = "1"

    logging.basicConfig(level=logging.INFO)

    args = parse_args()

    if args.subcommand == 'build':
        build(args)


if __name__ == "__main__":
    main()
