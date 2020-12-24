import pytest
import os
import filecmp
import shilofue.Doc as Doc
from shutil import rmtree, copyfile


# test files are in this directory
# future: change name in all other files as this one
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
    source_dir = os.path.join(test_source_dir, 'doc')
    mkdocs_dir = os.path.join(test_dir, 'doc')
    docs_dir = os.path.join(mkdocs_dir, 'docs')
    mkdocs_case_dir = os.path.join(docs_dir, 'foo')  # this should be generate by the code
    if os.path.isdir(mkdocs_dir):
        rmtree(mkdocs_dir)
    os.mkdir(mkdocs_dir)
    os.mkdir(docs_dir)
    os.mkdir(mkdocs_case_dir)

    # files that will be used in this test
    mkdocs_file = os.path.join(mkdocs_dir, 'mkdocs.yml')
    copyfile(os.path.join(source_dir, 'mkdocs.yml'), mkdocs_file)

    # call __init__ function
    myMkdoc = Doc.MKDOC(mkdocs_dir)

    #call __call__ function for a case
    # source directory for this case
    case_source_dir = os.path.join(source_dir, 'foo')
    myMkdoc('foo', os.path.join(case_source_dir), append_prm=True, images='DepthAverage')
    case_summary_file = os.path.join(mkdocs_case_dir, 'summary.md')
    case_prm_file = os.path.join(mkdocs_case_dir, 'case.prm')
    case_img_file = os.path.join(mkdocs_case_dir, 'img', 'DepthAverage_t0.00000000e+00.png')
    assert(os.path.isfile(case_summary_file))  # assert summary.md generated
    assert(os.path.isfile(case_img_file))  # assert hard link for image generated
    # assert file has the right content
    case_summary_file_std = os.path.join(case_source_dir, 'summary_std.md')
    assert(filecmp.cmp(case_summary_file, case_summary_file_std))
    
    # call __call__ function for a group
    group_source_dir = os.path.join(source_dir, 'foo_group')
    mkdocs_group_dir = os.path.join(docs_dir, 'foo_group')  # this should be generate by the code
    myMkdoc('foo_group', group_source_dir, append_prm=True, type='group', case_names=['foo1', 'foo2'], images='DepthAverage')
    group_summary_file = os.path.join(mkdocs_group_dir, 'summary.md')
    group_summary_file_std = os.path.join(group_source_dir, 'summary_std.md')
    assert(filecmp.cmp(group_summary_file, group_summary_file_std))
    # check for foo1
    case_source_dir = os.path.join(group_source_dir, 'foo1')
    mkdocs_case_dir = os.path.join(mkdocs_group_dir, 'foo1')
    case_summary_file = os.path.join(mkdocs_case_dir, 'summary.md')
    case_prm_file = os.path.join(mkdocs_case_dir, 'case.prm')
    # assert summary.md generated
    assert(os.path.isfile(case_summary_file))
    # assert file has the right content
    case_summary_file_std = os.path.join(case_source_dir, 'summary_std.md')
    assert(filecmp.cmp(case_summary_file, case_summary_file_std))
    # assert case.prm generated
    assert(os.path.isfile(case_prm_file))
    # check for foo2
    case_source_dir = os.path.join(group_source_dir, 'foo2')
    mkdocs_case_dir = os.path.join(mkdocs_group_dir, 'foo2')
    case_summary_file = os.path.join(mkdocs_case_dir, 'summary.md')
    case_prm_file = os.path.join(mkdocs_case_dir, 'case.prm')
    # assert summary.md generated
    assert(os.path.isfile(case_summary_file))
    # assert file has the right content
    case_summary_file_std = os.path.join(case_source_dir, 'summary_std.md')
    assert(filecmp.cmp(case_summary_file, case_summary_file_std))
    # assert case.prm generated
    assert(os.path.isfile(case_prm_file))

    # call __call__ function for a analysis
    case_dirs = ['foo', 'foo_group/foo1']
    images = ['DepthAverage']
    analysis_source_dir = os.path.join(source_dir, 'test_analysis')
    mkdocs_analysis_dir = os.path.join(docs_dir, 'test_analysis')  # this should be generate by the code
    # append an extra analysis of machine time
    # and Newton solver
    extra_analysis={'machine_time': {'step':2}, 'newton_solver': {}}
    myMkdoc('test_analysis', source_dir, append_prm=True, update=True, type='analysis', extra_analysis=extra_analysis,
            case_dirs=case_dirs, images=images)
    # myMkdoc('test_analysis', source_dir, append_prm=True, update=True, type='analysis', extra_analysis='machine_time',
    #        analyze_machine_step=2, case_dirs=case_dirs, images=images)
    # assertions
    analysis_summary_file = os.path.join(mkdocs_analysis_dir, 'summary.md')
    # assert summary.md generated
    assert(os.path.isfile(analysis_summary_file))
    analysis_summary_file_std = os.path.join(analysis_source_dir, 'summary_std.md')
    # assert file has the right content
    assert(filecmp.cmp(analysis_summary_file, analysis_summary_file_std))
    # assert images are linked
    analysis_image_dir = os.path.join(mkdocs_analysis_dir, 'img')
    images_ = ['foo_DepthAverage_t0.00000000e+00.png', 'foo_group-foo1_DepthAverage_t0.00000000e+00.png']
    for image_ in images_:
        assert(os.path.isfile(os.path.join(analysis_image_dir, image_)))
    # assert extra analysis is done
    assert(os.path.isfile(os.path.join(analysis_image_dir, "MachineTimeAnalysis.png")))
    assert(os.path.isfile(os.path.join(analysis_image_dir, "NewtonSolverAnalysis.png")))

    # assert mkdocs.yml file has the right contents
    mkdocs_file = os.path.join(mkdocs_dir, 'mkdocs.yml')
    mkdocs_file_std = os.path.join(source_dir, 'mkdocs_std.yml')
    assert(filecmp.cmp(mkdocs_file, mkdocs_file_std))

    # call __call__ function to update a new case in the previous group foo_group
    # todo
    group_source_dir = os.path.join(source_dir, 'foo_group')
    mkdocs_group_dir = os.path.join(docs_dir, 'foo_group')  # this should be generate by the code
    myMkdoc('foo_group', group_source_dir, append_prm=True, type='group', case_names=['foo1', 'foo2', 'foo3'], images='DepthAverage', update=True)
    # assert the summary file
    group_summary_file = os.path.join(mkdocs_group_dir, 'summary.md')
    group_summary_file_std = os.path.join(group_source_dir, 'summary_std.md')
    assert(filecmp.cmp(group_summary_file, group_summary_file_std))
    # assert the .yml file is updated
    mkdocs_file = os.path.join(mkdocs_dir, 'mkdocs.yml')
    mkdocs_file_std = os.path.join(source_dir, 'mkdocs_std1.yml')
    assert(filecmp.cmp(mkdocs_file, mkdocs_file_std))

            

