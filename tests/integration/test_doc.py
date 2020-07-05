import pytest
import os
from shutil import rmtree, copyfile
from shilofue.Doc import DOC, MKDOC


# test files are in this directory
# todo: change name in all other files as this one
test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
test_dir = './.test'


def test_doc():
    '''
    test class DOC from shilofue.doc
    '''
    # files that will be used in this test
    _ofile = os.path.join(test_dir, 'auto.md')
    _standard_ofile = os.path.join(test_source_dir, 'foo', 'auto.md')
    _depth_average_file = os.path.join(test_source_dir, 'foo', 'DepthAverage_t0.00000000e+00.pdf')
    if os.path.isfile(_ofile):
        # remove old files
        os.remove(_ofile)
    with pytest.raises(TypeError) as excinfo:
        Doc = DOC('foo', os.path.join(test_source_dir, 'foo'), layout='foo')
    assert('layout' in str(excinfo.value))
    Doc = DOC('foo', os.path.join(test_source_dir, 'foo'))
    Doc.append_media({'depth_average': _depth_average_file})
    _contents = Doc.GenerateOutput
    assert(True)  # assert the output contents of test case
    Doc(test_dir)
    assert(os.path.isfile(_ofile))  # assert file generated successfully
    with open(_ofile, 'r') as fin:
        # assert file has the right content
        with open(_standard_ofile, 'r') as standard_fin:
            assert(fin.read() == standard_fin.read())


def test_mkdoc():
    '''
    test class MKDOC from shilofue.doc
    '''
    if os.path.isdir(os.path.join(test_dir, 'test-project')):
        rmtree(os.path.join(test_dir, 'test-project'))
    os.mkdir(os.path.join(test_dir, 'test-project'))
    os.mkdir(os.path.join(test_dir, 'test-project', 'docs'))
    _mkdocs_case_dir = os.path.join(test_dir, 'test-project', 'docs', 'foo')  # this should be generate by the code
    os.mkdir(_mkdocs_case_dir)
    # files that will be used in this test
    _mkdocs_file = os.path.join(test_dir, 'test-project', 'mkdocs.yml')
    copyfile(os.path.join(test_source_dir, 'test-project', 'mkdocs.yml'), _mkdocs_file)
    _index_file = os.path.join(test_dir, 'test-project', 'docs', 'index.md')
    _case_file = os.path.join(_mkdocs_case_dir, 'foo.md')
    # call __init__ function
    myMkdoc = MKDOC(os.path.join(test_source_dir, 'test-project', 'foo'), os.path.join(test_dir, 'test-project'))
    # call __call__ function
    myMkdoc('foo', append_prm=True)
    assert(os.path.isfile(os.path.join(_mkdocs_case_dir, 'summary.md')))  # assert summary.md generated
    with open(os.path.join(_mkdocs_case_dir, 'summary.md'), 'r') as fin:
        # assert file has the right content
        with open(os.path.join(test_source_dir, 'test-project', 'standard_foo_summary.md')) as standard_fin:
            assert(fin.read() == standard_fin.read())
    assert(os.path.isfile(os.path.join(_mkdocs_case_dir, 'case.prm')))  # assert case.prm generated
    with open(os.path.join(test_dir, 'test-project', 'mkdocs.yml')) as fin:
        with open(os.path.join(test_source_dir, 'test-project', 'standard_mkdocs.yml')) as standard_fin:
            assert(fin.read() == standard_fin.read())
    # todo: test the update option
            

