import os
import json
import filecmp
import sys
import shilofue.TwoDSubduction as TwoDSubduction
import shilofue.Parse as Parse
import shilofue.ParsePrm as ParsePrm
from shutil import rmtree

ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
test_dir = '.test'
project_pp_json = os.path.join(ASPECT_LAB_DIR, 'files', 'TwoDSubduction', 'post_process.json')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


def test_bash_options():
    """
    test BASH_OPTIONS class
    """
    pass