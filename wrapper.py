import os
import sys
import argparse

def check_arg(args=None):
    """
    Function to parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="ADD TITLE OF SCRIPT HERE (shows on help)")
    parser.add_argument("-r", "--reference",
    help="reference transcriptome",
    required=True)
    parser.add_argument("-c", "--condensed",
    help="condensed pipeline", 
    required=True)
    return parser.parse_args(args)

def main():
    pass

if __name__ == "__main__":
    main()