import pytest
import os
import shilofue.Doc as Doc
from shutil import rmtree, copyfile


# test files are in this directory
# todo: change name in all other files as this one
test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
test_dir = './.test'


def test_doc():
    '''
    test class DDOC from shilofue.Doc
    '''
    # files that will be used in this test
    _ofile = os.path.join(test_dir, 'auto.md')
    _standard_ofile = os.path.join(test_source_dir, 'foo', 'auto.md')
    _depth_average_file = os.path.join(test_source_dir, 'foo', 'DepthAverage_t0.00000000e+00.pdf')
    if os.path.isfile(_ofile):
        # remove old files
        os.remove(_ofile)
    with pytest.raises(TypeError) as excinfo:
        Ddoc = Doc.DDOC('foo', os.path.join(test_source_dir, 'foo'), layout='foo')
    assert('layout' in str(excinfo.value))
    Ddoc = Doc.DDOC('foo', os.path.join(test_source_dir, 'foo'))
    Ddoc.append_media({'depth_average': _depth_average_file})
    _contents = Ddoc.GenerateOutput
    assert(True)  # assert the output contents of test case
    Ddoc(test_dir)
    assert(os.path.isfile(_ofile))  # assert file generated successfully
    with open(_ofile, 'r') as fin:
        # assert file has the right content
        with open(_standard_ofile, 'r') as standard_fin:
            assert(fin.read() == standard_fin.read())


def test_mkdoc():
    '''
    test class MKDOC from shilofue.doc
    '''
    _mkdocs_dir = os.path.join(test_dir, 'test-project')
    _docs_dir = os.path.join(_mkdocs_dir, 'docs')
    _mkdocs_case_dir = os.path.join(_docs_dir, 'foo')  # this should be generate by the code
    if os.path.isdir(_mkdocs_dir):
        rmtree(_mkdocs_dir)
    os.mkdir(_mkdocs_dir)
    os.mkdir(_docs_dir)
    os.mkdir(_mkdocs_case_dir)
    # files that will be used in this test
    _mkdocs_file = os.path.join(_mkdocs_dir, 'mkdocs.yml')
    copyfile(os.path.join(test_source_dir, 'test-project', 'mkdocs.yml'), _mkdocs_file)
    _index_file = os.path.join(_docs_dir, 'index.md')
    _case_file = os.path.join(_mkdocs_case_dir, 'foo.md')
    # call __init__ function
    myMkdoc = Doc.MKDOC(_mkdocs_dir, images='DepthAverage')
    ############# call __call__ function for a case ###############
    myMkdoc('foo', os.path.join(test_source_dir, 'test-project', 'foo'), append_prm=True)
    _case_summary_file = os.path.join(_mkdocs_case_dir, 'summary.md')
    _case_prm_file = os.path.join(_mkdocs_case_dir, 'case.prm')
    _case_img_file = os.path.join(_mkdocs_case_dir, 'img', 'DepthAverage_t0.00000000e+00.pdf')
    assert(os.path.isfile(_case_summary_file))  # assert summary.md generated
    assert(os.path.isfile(_case_img_file))  # assert hard link for image generated
    # assert file has the right content
    with open(_case_summary_file, 'r') as fin:
        _case_summary_contents = fin.read()
    with open(os.path.join(test_source_dir, 'test-project', 'standard_foo_summary.md')) as standard_fin:
        assert(_case_summary_contents == standard_fin.read())
    assert(os.path.isfile(_case_prm_file))  # assert case.prm generated
    ############# call __call__ function for a group ##############
    _mkdocs_group_dir = os.path.join(_docs_dir, 'foo_group')  # this should be generate by the code
    myMkdoc('foo_group', os.path.join(test_source_dir, 'test-project', 'foo_group'), append_prm=True, type='group', case_names=['foo1', 'foo2'])
    _group_summary_file = os.path.join(_mkdocs_group_dir, 'summary.md')
    # check group summary, todo
    # check for foo1
    _mkdocs_subcase_dir = os.path.join(_mkdocs_group_dir, 'foo1')
    _case_summary_file = os.path.join(_mkdocs_subcase_dir, 'summary.md')
    _case_prm_file = os.path.join(_mkdocs_subcase_dir, 'case.prm')
    assert(os.path.isfile(_case_summary_file))  # assert summary.md generated
    # assert file has the right content
    with open(_case_summary_file, 'r') as fin:
        _case_summary_contents = fin.read()
    with open(os.path.join(test_source_dir, 'test-project', 'standard_foo1_summary.md')) as standard_fin:
        assert(_case_summary_contents == standard_fin.read())
    assert(os.path.isfile(_case_prm_file))  # assert case.prm generated
    # check for foo2
    _mkdocs_subcase_dir = os.path.join(_mkdocs_group_dir, 'foo2')
    _case_summary_file = os.path.join(_mkdocs_subcase_dir, 'summary.md')
    _case_prm_file = os.path.join(_mkdocs_subcase_dir, 'case.prm')
    assert(os.path.isfile(_case_summary_file))  # assert summary.md generated
    # assert file has the right content
    with open(_case_summary_file, 'r') as fin:
        _case_summary_contents = fin.read()
    with open(os.path.join(test_source_dir, 'test-project', 'standard_foo2_summary.md')) as standard_fin:
        assert(_case_summary_contents == standard_fin.read())
    assert(os.path.isfile(_case_prm_file))  # assert case.prm generated
    ############ assert mkdocs.yml file has the right contents ########3
    with open(os.path.join(test_dir, 'test-project', 'mkdocs.yml')) as fin:
        _yml_contents = fin.read()
    with open(os.path.join(test_source_dir, 'test-project', 'standard_mkdocs.yml')) as standard_fin:
        assert(_yml_contents == standard_fin.read())
    # todo_future: test the update option
            

