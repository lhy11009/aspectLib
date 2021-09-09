import sys, os
import argparse
import shilofue

def main():
    """
    main function
    """
    # parse options and arguments
    parser = argparse.ArgumentParser(description='todo')
    parser.add_argument('-i', '--inputs', type=str,
                        default = None,
                        help='some inputs')
    # parse
    arg = parser.parse_args()

    # todo
    pass


if __name__ == "__main__":
    main()