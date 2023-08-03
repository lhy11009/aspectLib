import numpy
import json
import os, sys
import re
import filecmp
from shilofue.CaseOptions import ASPECT_LAB_DIR
import shilofue.json_files
import pdb
import shilofue.Plot as Plot
from matplotlib import pyplot as plt
from matplotlib import cm
from shutil import copyfile
from pathlib import Path

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
JSON_FILE_DIR = os.path.join(ASPECT_LAB_DIR, "files", "json_examples")

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
from Utilities import my_assert, re_neat_word, re_count_indent, touch



class DDOC():
    '''
    class for documentation
    Attributes:
        case_name(str) - case name
        dir(str) - directory of case data
        layout(dict) - layout of doc
        media(dict) - media name: route
    '''
    def __init__(self, _case_name, _dir, **kwargs):
        '''
        Inputs:
            _case_name(str) - case name
            _dir(str) - directory of case data
            _kwargs:
                layout(dict):
                    layout of doc
                layout_json(str):
                    file name to read layout from
        '''
        self.case_name = _case_name
        self.idir = _dir
        self.media = {}
        # get layout
        try:
            _layout = kwargs['layout']
            my_assert(type(_layout) == dict, TypeError, "layout must be a dictionary")
            self.layout = _layout
        except KeyError:
            _layout_json = kwargs.get('layout_json', 'DocLayout.json')
            json_path = os.path.join(JSON_FILE_DIR, _layout_json)
            with open(json_path, 'r') as fin:
                _layout = json.load(fin)
            my_assert(type(_layout) == dict, TypeError, "layout must be a dictionary")
            self.layout = _layout

    def append_media(self, _media):
        '''
        Inputs: 
            _media(dict) - name of media: route
        '''
        for key, value in _media.items():
            self.media[key] = value

    def __call__(self, _odir, **kwargs):
        '''
        call and generate output md file
        Inputs:
            _odir(str) - route of output
            kwargs:
                filename - name of output file
        '''
        _ofile = kwargs.get('filename', 'auto.md')
        if not os.path.isdir(_odir):
            os.mkdir(_odir)
        _ofile = os.path.join(_odir, _ofile)
        _contents = self.GenerateOutput()
        with open(_ofile, 'w') as fout:
            fout.write(_contents)

    def GenerateOutput(self, **kwargs):
        '''
        Generate contents of output file
        Outputs:
            contents(str): output contents
        '''
        _contents = ''
        for key, value in self.layout.items():
            if key == 'case_name':
                _contents += '# %s\n' % self.case_name
            else:
                _contents += '# %s\n\n' % key
                if type(value) == str:
                    try:
                        _contents += (ConvertMediaMKD(value, self.media[value]) + '\n')
                    except KeyError:
                        _contents += 'None\n'
                elif type(value) == list:
                    for _sub_value in value:
                        my_assert(type(_sub_value) == str, TypeError, "value must be string or a list of string")
                        _contents += "## %s \n\n" % _sub_value
                        try:
                            _contents += (ConvertMediaMKD(_sub_value, self.media[_sub_value]) + '\n')
                        except KeyError:
                            _contents += 'None\n\n'
            _contents += '\n'
        return _contents
            

def ConvertMediaMKD(_name, _filename):
    '''
    return a media string in markdown
    Inputs:
        _name(str): name of media
        _filename(str): filename of media
    Outputs:
        _output(str): media string
    '''
    _output = "![%s](%s)" % (_name, _filename)
    return _output


def ConvertLinkMKD(_name, _filename):
    '''
    return a link string in markdown
    Inputs:
        _name(str): name of media
        _filename(str): filename of media
    Outputs:
        _output(str): media string
    '''
    _output = "[%s](%s)" % (_name, _filename)
    return _output


class MKDOC():
    '''
    A class to use with mkdocs
    Attributes:
        name(str) - case name
        dir(str) - directory of case data
        imgs(list) - list of images
        new_files(dict) - name: new file
    '''
    def __init__(self, _odir, **kwargs):
        '''
        Inputs:
            _odir(str) - directory of mkdocs project
        '''
        self.odir = _odir
        self.new_files = {}

    def AttachImage(self, _img):
        '''
        Inputs:
            _img(list or str) - img files
        '''
        if type(_img) == str:
            self.imgs.append(_img)
        elif type(_img) == list:
            self.imgs += _img
        else:
            raise TypeError('img must by a str or a list')

    def __call__(self, _name, _dir, **kwargs):
        '''
        Call function and write to file that will be used bu mkdoc
        Inputs:
            kwargs: 
                update(bool) - if update an existing case
                append_prm(bool) - if append a prm file
                prm(str) - name of a prm file
                extra_images(list) - names of extra images to append
        '''
        self.new_files = {}
        self.imgs = kwargs.get('images', [])
        update = kwargs.get('update', False)
        append_prm = kwargs.get('append_prm', False)
        _prm = kwargs.get('prm', 'case.prm')
        _type = kwargs.get('type', 'case')

        # These directory and files need to be preexist
        assert(os.path.isdir(self.odir))
        assert(os.path.isdir(os.path.join(self.odir, 'docs')))
        _mcdocs_file = os.path.join(self.odir, 'mkdocs.yml')
        assert(os.path.isfile(_mcdocs_file))
        _index_file = os.path.join(self.odir, 'docs', 'index.md')

        # deal with different types of entry
        _target_dir = os.path.join(self.odir, 'docs', _name)
        if _type == 'case':
            # append a case
            self.AppendCase(_name, _dir, _target_dir, update=update, append_prm=append_prm)
        elif _type == 'group':
            # append a group
            _case_names = kwargs.get('case_names', None)
            my_assert(_case_names is not None, ValueError, 'For a group, case names cannot be None. Valid names must be given')
            self.AppendGroup(_name, _dir, _case_names, _target_dir, update=update, append_prm=append_prm)
        elif _type == 'analysis':
            # append an analysis
            _case_dirs = kwargs.get('case_dirs', [])
            # here the _dir is the project directory and case_dirs are relative directories to that
            self.AppendAnalysis(_name, _dir, _case_dirs, _target_dir, kwargs)
        else:
            raise ValueError("Type must be 'case', 'group', 'analysis' ")
        if self.new_files != {}:
            self.RenewMkdocsYml(_name)
        
    def AppendCase(self, _name, _dir, _target_dir, **kwargs):
        '''
        Append a case to doc
        Inputs:
            _target_dir(str): directory to put outputs
            kwargs:
                update(bool) - if update an existing case
                append_prm(bool) - if append a prm file
                prm(str) - name of a prm file
                basename(str) - base name that goes beform path of the case in .yml file,
                                used to figure out a case of a group
        '''
        update = kwargs.get('update', False)
        append_prm = kwargs.get('append_prm', False)
        _prm = kwargs.get('prm', 'case.prm')
        _base_name = kwargs.get('basename', None)
        if not os.path.isdir(_target_dir):
            os.mkdir(_target_dir)
        # get images
        _img_dir = os.path.join(_dir, 'img')
        _imgs = ReturnFileList(_img_dir, self.imgs)
        # create hard links in target_dir
        _target_img_dir = os.path.join(_target_dir, 'img')
        if not os.path.isdir(_target_img_dir):
            os.mkdir(_target_img_dir)
        for _img in _imgs:
            _file = os.path.join(_img_dir, _img)
            _target_file = os.path.join(_target_img_dir, _img)
            if not os.path.isfile(_target_file):
                os.link(_file, _target_file)
            elif filecmp.cmp(_file, _target_file) is False:
                os.remove(_target_file)
                os.link(_file, _target_file)
        # Append a summary.md
        if os.path.isfile(os.path.join(_target_dir, 'summary.md')):
            if update == True:
                self.GenerateCaseMkd(_dir, _target_dir, update=True, images=_imgs)
        else:
            _filename = self.GenerateCaseMkd(_dir, _target_dir, images=_imgs)
            # in a mkdocs file, files are listed as 'name/_filename'
            _summary = os.path.join(_name, os.path.basename(_filename))
            if _base_name is None:
                self.new_files['Summary'] = _summary
            else:
                # a subcase of a group
                try:
                    self.new_files[_name]['Summary'] = os.path.join(_base_name, _summary)
                except KeyError:
                    self.new_files[_name] = {}
                    self.new_files[_name]['Summary'] = os.path.join(_base_name, _summary)
        # Append a prm file if the append_prm option is True
        if append_prm:
            if os.path.isfile(os.path.join(_target_dir, _prm)):
                if update == True:
                    _prm_source_file = os.path.join(_dir, _prm)
                    assert(os.path.isfile(_prm_source_file))
                    copyfile(_prm_source_file, os.path.join(_target_dir, _prm))
            else:
                _prm_source_file = os.path.join(_dir, _prm)
                assert(os.path.isfile(_prm_source_file))
                copyfile(_prm_source_file, os.path.join(_target_dir, _prm))
                _parameters = os.path.join(_name, _prm)
                if _base_name is None:
                    self.new_files['Parameters'] = _parameters
                else:
                    # a subcase of a group
                    try:
                        self.new_files[_name]['Parameters'] = os.path.join(_base_name, _parameters)
                    except KeyError:
                        self.new_files[_name] = {}
                        self.new_files[_name]['Parameters'] = os.path.join(_base_name, _parameters)
    
    def AppendGroup(self, _name, _dir, _case_names, _target_dir, **kwargs):
        '''
        Append a group to doc
        Inputs:
            _case_names(list): names of cases within this group
            _target_dir(str): directory to put outputs
            kwargs:
                update(bool) - if update an existing case
                append_prm(bool) - if append a prm file
                prm(str) - name of a prm file
        '''
        update = kwargs.get('update', False)
        append_prm = kwargs.get('append_prm', False)
        _prm = kwargs.get('prm', 'case.prm')
        if not os.path.isdir(_target_dir):
            os.mkdir(_target_dir)
        # Append a summary.md
        if os.path.isfile(os.path.join(_target_dir, 'summary.md')):
            if update == True:
                self.GenerateGroupMkd(_dir, _target_dir, update=True)
        else:
            _filename = self.GenerateGroupMkd(_dir, _target_dir)
            # in a mkdocs file, files are listed as 'name/_filename'
            self.new_files['Summary'] = os.path.join(_name, os.path.basename(_filename))
        for _case_name in _case_names:
            _case_dir = os.path.join(_dir, _case_name)
            _case_target_dir = os.path.join(_target_dir, _case_name)
            self.AppendCase(_case_name, _case_dir, _case_target_dir, update=update, append_prm=append_prm, basename=_name)

    def AnalyzeExtra(self, _name, _project_dir, _case_dirs, _target_dir, extra_analysis, kwargs):
        # extra procedures in an analysis
        # Inputs:
        #   extra_analysis(str): type of extra analysis
        #   kwargs(dict): dictionary of options
        for key,value in extra_analysis.items():
            if key == 'machine_time':
                # todo change to a class 
                my_assert(type(value)==dict, TypeError, "AnalyzeExtra: settings to an option(a key in the analysis dict) must be a dict")
                AnalyzeMachineTime = ANALYZEMACHINETIME(_project_dir, _case_dirs, _target_dir, value)
                AnalyzeMachineTime(value)
            if key == 'newton_solver':
                my_assert(type(value)==dict, TypeError, "AnalyzeExtra: settings to an option(a key in the analysis dict) must be a dict")
                self.AnalyzeNewtonSolver(_project_dir, _case_dirs, _target_dir, value)

    
    def AnalyzeNewtonSolver(self, _project_dir, _case_dirs, _target_dir, kwargs):
        '''
        generate a comparison of solver output
        '''
        # Instantiate an object for Solver output
        NewtonSolver = Plot.NEWTON_SOLVER_PLOT('NewtonSolver')
        
        # Initialize plot
        fig, ax = plt.subplots()

        # create a color table
        normalizer = [float(i)/(len(_case_dirs)-1) for i in range(len(_case_dirs))]
        colors = cm.rainbow(normalizer)

        # loop for cases, read data and plot
        for i in range(len(_case_dirs)):
            _dir = _case_dirs[i]
            _img_dir = os.path.join(_project_dir, _dir, 'img')
            _output_dir = os.path.join(_project_dir, _dir, 'output')
            solver_output_file = os.path.join(_output_dir, 'solver_output')
            if os.path.isfile(solver_output_file):
                NewtonSolver.ReadHeader(solver_output_file)  # inteprate header information
        
                # catch possible vacant file
                state = NewtonSolver.ReadData(solver_output_file)  # read data
                if state == 1:
                    continue
        
                # manage output data for all steps
                _data_list = NewtonSolver.ManageDataAll()

                # plot
                options={'xname': 'Number_of_nonlinear_iteration', 'yname': 'Relative_nonlinear_residual', 'log_y': 1,
                        'line': '.-', 'color': colors[i], 'label': _dir}
                NewtonSolver.Plot(_data_list, ax, options)
        
        # save figure
        _target_img_dir = os.path.join(_target_dir, 'img')
        fileout = os.path.join(_target_img_dir, 'NewtonSolverAnalysis.png')
        ax.set(xlabel="Number of nonlinear iteration", ylabel="Relative nonlinear residual")
        ax.legend()
        ax.set_title("Relative Nonlinear Residual")
        fig.savefig(fileout)
    
    def AppendAnalysis(self, _name, _project_dir, _case_dirs, _target_dir, kwargs):
        '''
        Append a analysis to doc
        Inputs:
            _name(str): name of analysis
            _case_dirs(list): list of directory of cases to include
            _target_dir(str): directory to put outputs
        '''
        update = kwargs.get('update', False)

        # check on target directory 
        if not os.path.isdir(_target_dir):
            os.mkdir(_target_dir)
        
        # create hard link for images
        _target_img_dir = os.path.join(_target_dir, 'img')
        if not os.path.isdir(_target_img_dir):
            os.mkdir(_target_img_dir)
        # loop imgs first, so as to create a two-d list
        _imgs_list = []
        for i in range(len(self.imgs)):
            _imgs_list.append([])
        # _imgs_list = [[]]*len(self.imgs)
        for i in range(len(self.imgs)):
            img = self.imgs[i]
            for _dir in _case_dirs:
                _img_dir = os.path.join(_project_dir, _dir, 'img')
                _imgs = ReturnFileList(_img_dir, [img])
                # transfer _dir to a name to append
                _dir_transfered = re.sub(os.sep, '-', _dir)
                #_dir_transfered = os.path.basename(_dir)
                # create hard links in target_dir
                for _img in _imgs:
                    _file = os.path.join(_img_dir, _img)
                    _target_file = os.path.join(_target_img_dir, "%s_%s" %(_dir_transfered, _img))
                    if not os.path.isfile(_target_file):
                        os.link(_file, _target_file)
                    elif filecmp.cmp(_file, _target_file) is False:
                        os.remove(_target_file)
                        os.link(_file, _target_file)
                    _imgs_list[i].append(_target_file)

        # deal with extra analysis
        extra_analysis = kwargs.get('extra_analysis', {})
        my_assert(type(extra_analysis)==dict, TypeError, "AppendAnalysis: extra_analysis must be a dict")
        self.AnalyzeExtra(_name, _project_dir, _case_dirs, _target_dir, extra_analysis, kwargs)
        
        # Append a summary.md
        # append image information
        _base_name = kwargs.get('basename', None)
        if os.path.isfile(os.path.join(_target_dir, 'summary.md')):
            if update == True:
                _filename = self.GenerateAnalysisMkd(_target_dir, _case_dirs, images=_imgs_list, extra_analysis=extra_analysis)
        else:
            _filename = self.GenerateAnalysisMkd(_target_dir, _case_dirs, images=_imgs_list, extra_analysis=extra_analysis)
            # in a mkdocs file, files are listed as 'name/_filename'
            _summary = os.path.join(_name, os.path.basename(_filename))
            if _base_name is None:
                self.new_files['Summary'] = _summary
            else:
                # a subcase of a group
                self.new_files[_name]['Summary'] = os.path.join(_base_name, _summary)
        pass

    def GenerateCaseMkd(self, _dir, _target_dir, **kwargs):
        '''
        Generate markdown file of a case
        Inputs:
            _dir(str): directory of case
            _target_dir(str): directory of this case
            kwargs:
                filename(str): name of the file
                images(list of str): images to append
        Returns:
            _filename(str): file generated
        '''
        _filename = kwargs.get('filename', 'summary.md')
        _imgs = kwargs.get('images', [])
        _filename = os.path.join(_target_dir, _filename)
        _auto_mkd_file = os.path.join(_dir, 'auto.md')
        _extra_mkd_file = os.path.join(_dir, 'extra.md')
        
        # copy the auto.mkd file from case directory
        # my_assert (os.path.isfile(_auto_mkd_file), AssertionError, "%s is unreadable" % _auto_mkd_file)
        # copyfile(_auto_mkd_file, _filename)
        
        if os.path.isfile(_auto_mkd_file):
            copyfile(_auto_mkd_file, _filename)
        else:
            touch(_filename)
            
        
        # apppend contents of extra.mkd at the end of case.mkd
        if (os.path.isfile(_extra_mkd_file)):
            with open(_extra_mkd_file, 'r') as fin:
                _contents = fin.read()
            with open(_filename, 'a') as fout:
                fout.write('\n')
                fout.write(_contents)
        
        # append images
        if _imgs != []:
            _contents = self.GenerateImageMkd(_dir, _imgs)
            with open(_filename, 'a') as fout:
                fout.write('\n')
                fout.write(_contents)
        return _filename

    def GenerateImageMkd(self, _dir, _files, **kwargs):
        '''
        Generate markdown file content for figures
        Inputs:
            _dir(str): directory of case
            _files(list of str): files to be included
        Return:
            _contents: contents of markdown file
        '''
        _contents = '# Plots\n\n'
        for _file in _files:
            _contents += '## %s\n\n' % _file
            _relative_route = os.path.join('img', _file)
            _line = '![%s](%s)' %  (_file, _relative_route)
            _contents += '%s\n\n' % (_line)
        return _contents
    
    def GenerateGroupMkd(self, _dir, _target_dir, **kwargs):
        '''
        Generate markdown file of a group
        Inputs:
            _dir(str): directory of case
            _target_dir(str): directory of this case
            kwargs:
                filename(str): name of the file
        Returns:
            _filename(str): file generated
        '''
        _filename = kwargs.get('filename', 'summary.md')
        _filename = os.path.join(_target_dir, _filename)
        _auto_mkd_file = os.path.join(_dir, 'auto.md')
        _extra_mkd_file = os.path.join(_dir, 'extra.md')
        assert (os.path.isfile(_auto_mkd_file))
        copyfile(_auto_mkd_file, _filename)  # copy the auto.mkd file from case directory
        if (os.path.isfile(_extra_mkd_file)):
            # apppend contents of extra.mkd at the end of case.mkd
            with open(_extra_mkd_file, 'r') as fin:
                _contents = fin.read()
            with open(_filename, 'a') as fout:
                fout.write(_contents)
        return _filename
                
        
    def GenerateAnalysisMkd(self, _target_dir, _case_dirs, **kwargs):
        '''
        Generate markdown file of a analysis
        Inputs:
            _target_dir(str): directory of this case
            kwargs:
                filename(str): name of the file
        Returns:
            _filename(str): file generated
        '''
        # get filename
        _filename = kwargs.get('filename', 'summary.md')
        _filename = os.path.join(_target_dir, _filename)

        # write file header
        contents = ''
        _name = os.path.basename(_target_dir)
        contents += '# Analysis %s\n\n' % _name
        
        # append names of cases
        contents += '## This includes cases: \n\n'
        for _case_dir in _case_dirs:
            _case_link = os.path.join('..', _case_dir, 'summary.md')
            contents += '* %s\n' % ConvertLinkMKD(_case_dir, _case_link)
        contents += '\n'
       
        # append imgs
        imgs_ = kwargs.get('images', [[]]*len(self.imgs))
        my_assert(len(self.imgs)==len(imgs_), ValueError,
                  "GenerateAnalysisMkd: images input must be a list that have the same length with self.imgs")
        for i in range(len(self.imgs)):
            img_root = self.imgs[i]
            contents += '## Image %s\n\n' % img_root
            for img_ in imgs_[i]:
                basename_ = os.path.basename(img_)
                relative_route = os.path.join('img', os.path.basename(basename_))
                contents += '%s\n' % basename_
                contents += '%s\n\n' % ConvertMediaMKD(basename_, relative_route)

        # append extra analysis
        extra_analysis = kwargs.get('extra_analysis', {})
        for key, value in extra_analysis.items():
            contents += '## %s\n\n' % key
            if key == 'machine_time':
                contents += 'Here we show machine time (core hrs) for each case\n'
                contents += '%s\n\n' % ConvertMediaMKD("MachineTimeAnalysis.png", "img/MachineTimeAnalysis.png")
            if key == 'newton_solver':
                contents += 'Here we show solver output for each case\n'
                contents += '%s\n\n' % ConvertMediaMKD("NewtonSolverAnalysis.png", "img/NewtonSolverAnalysis.png")

        # write
        with open(_filename, 'w') as fout:
            fout.write(contents)
        pass

        return _filename
        

    def RenewMkdocsYml(self, _name):
        '''
        Renew the mkdocs.yml file in the project directory
        '''
        _filename = os.path.join(self.odir, 'mkdocs.yml')
        _start = None
        _end = None
        with open(_filename, 'r') as fin:
            # read in the old file
            _lines = fin.read().split('\n')
        for i in range(len(_lines)):
            # read in the nav part
            _line = _lines[i]
            if _start is None and re.match('^nav', _line):
                _start = i + 1
                _previous_indent = re_count_indent(_lines[i])
                _indent = re_count_indent(_lines[i+1])
            elif _start is not None and re_count_indent(_line) == _previous_indent:
                _end = i
                break
        my_assert(_start is not None and _end is not None, TypeError, 'Cannot find start and end of the nav part')
        _nav_dict, _temp= ExtractNav(_lines[_start: _end], previous=_indent)
        assert(_temp == _end - _start - 1)
        try:
            # this case is already in the nav part of the yml file
            value = _nav_dict[_name]
            my_assert(type(value) == dict, TypeError,
                      'entry for a case must be a single dict, prepared to include dictionary in the future')
            _nav_dict[_name] = {**value, **self.new_files}  # merge and substitute value in first dict with value in second dict
        except KeyError:
            # this case is new to the yml file
            _nav_dict[_name] = self.new_files
        _new_lines = _lines[0: _start]
        _new_lines += ProduceNav(_nav_dict)
        _new_lines += _lines[_end: -1]
        with open(_filename, 'w') as fout:
            for _line in _new_lines:
                fout.write(_line + '\n')


class ANALYZEMACHINETIME():
    '''
    Analyze machine time outputs
    Atributes:
        project_dir(str) : directory of project
        case_dirs(str) : directory of cases included
        target_dir: directory of mkdocs where the results are kept
    ''' 
    def __init__(self, project_dir, case_dirs, target_dir, kwargs={}):
        '''
        Inputs:
            kwargs(dict): options, optional
        '''
        self.project_dir = project_dir
        self.case_dirs = case_dirs
        self.target_dir = target_dir
        # Instantiate an object for MachineTime
        self.MachineTime = Plot.MACHINE_TIME_PLOT('MachineTime')

    def __call__(self, kwargs):
        '''
        Analysis and plot
        '''
        my_assert(type(kwargs) == dict, TypeError, "ANALYZEMACHIENTIME: __call__: kwargs must be a dictionary")

        # read data
        step = kwargs.get('step', 1)
        self.ReadData(step)

        # do analysis
        type_ = kwargs.get('type', 'strong_scaling')
        if type_ == 'strong_scaling':
            self.StrongScaling(step, kwargs)
        else:
            raise ValueError("type_ could only be one of: 'strong_scaling'")
    
    def ReadData(self, step):
        '''
        import data for analysis
        '''
        self.machine_times = []
        self.number_of_cpus = []
        for _dir in self.case_dirs:
            _img_dir = os.path.join(self.project_dir, _dir, 'img')
            _output_dir = os.path.join(self.project_dir, _dir, 'output')
            machine_time_file = os.path.join(_output_dir, 'machine_time')
            if os.path.isfile(machine_time_file):
                machine_time, number_of_cpu = self.MachineTime.GetStepMT(machine_time_file, step)
                self.machine_times.append(machine_time)
                self.number_of_cpus.append(number_of_cpu)
            else:
                self.machine_times.append(None)
                self.number_of_cpus.append(None)

    def StrongScaling(self, step, kwargs={}):
        '''
        todo
        generate a strong scaling plot of wall-time vs cores
        Inputs:
            kwargs(dict): options
        '''
        # plot
        _target_img_dir = os.path.join(self.target_dir, 'img')
        fileout_png = os.path.join(_target_img_dir, 'MachineTimeAnalysis.png')
        fileout_pdf = os.path.join(_target_img_dir, 'MachineTimeAnalysis.pdf')
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))

        # create a color table
        normalizer = [float(i)/(len(self.case_dirs)-1) for i in range(len(self.case_dirs))]
        colors = cm.rainbow(normalizer)

        # generate plot
        for i in range(len(self.machine_times)):
            if self.number_of_cpus[i] != None and self.machine_times[i] != None:
                axs[0].loglog(self.number_of_cpus[i], self.machine_times[i]*3600.0, 'o', color=colors[i], label=self.case_dirs[i])
                axs[1].loglog(self.number_of_cpus[i], self.machine_times[i]/self.number_of_cpus[i]*3600.0, 'o', color=colors[i], label=self.case_dirs[i])
        axs[0].set(xlabel="Number of CPU", ylabel="Machine Time (core*second)")
        axs[1].set(xlabel="Number of CPU", ylabel="Wall-clock Time (second)")
        axs[0].legend()
        axs[0].set_title("Machine Time at step %d" % step)
        axs[1].set_title("Wall-clock Time at step %d" % step)
        fig.tight_layout()
        fig.savefig(fileout_png)
        # attach a pdf
        attach_pdf = kwargs.get('attach_pdf', 0)
        if attach_pdf:
            fig.savefig(fileout_pdf)


def ExtractNav(_lines, **kwargs):
    '''
    extract the nav message from plain text
    Inputs:
        _contents(str): plain text contents
        kwargs:
        fix the name with the produce function
    '''
    _previous = kwargs.get('previous', 0)
    _at = kwargs.get('at', 0)
    assert(_at <= len(_lines))
    _odict = {}
    i = -1
    for i in range(_at, len(_lines)):
        _line = _lines[i]
        my_assert('-' in _line, TypeError, "Each line must start with \'-\'")
        _now = re_count_indent(_line)
        if _now < _previous:
            # this level is all done
            break
        elif _now > _previous:
            # come to lower level, which is already fixed by interation
            continue
        else:
            # still on this level, add new memember to dictionary
            key, value = SeparateNavPattern(_line)
            if value == '':
                # going to the next level, call iteration
                _next = re_count_indent(_lines[i+1])
                my_assert(_next > _now, ValueError, 'A group with no cases or a Case with no files are not supported')
                _odict[key], i = ExtractNav(_lines, previous=_next, at=i+1)
            else:
                # add key and value 
                _odict[key] = value
    return _odict, i


def ProduceNav(_idict, **kwargs):
    '''
    Produce the nav part of the mkdocs from a dictionary
    Inputs:
        _idict(dict): dictionary that has the information of nav
        kwargs:
            now(int): how many indent we have now. This is used to do iteration
            indent(int): indent between two adjacent levels
    Outputs:
        _lines(list of str): contents of nav in the .yml file
    '''
    _now = kwargs.get('now', 4)
    _indent = kwargs.get('indent', 4)  # indent between adjacent levels
    _lines = []
    for key, value in _idict.items():
        _line = ' '*_now + '- ' + key + ':'
        if type(value) == str:
            _line += (' ' + value)
            _lines.append(_line)  # append this line to the list
        elif type(value) == dict:
             _lines.append(_line)  # append this line to the list
             _lines += ProduceNav(value, now=_now+_indent, indent=_indent)
        else:
             raise TypeError('Value must be str or dict')
    return _lines


def SeparateNavPattern(_pattern):
    '''
    separet a Nav patttern and return key and value
        Inputs:
            _pattern(str): a nav patter, "- key: value"
    '''
    _temp = _pattern.split('-', maxsplit=1)[1]
    key = re_neat_word(_temp.split(':')[0])
    value = re_neat_word(_temp.split(':')[1])
    return key, value


def UpdateProjectDoc(_project_dict, _project_dir, **kwargs):
    '''
    Update doc for all cases in this project
    Inputs:
        kwargs(dict): options
            analysis: a dictionary of tests
    '''
    _mkdocs = kwargs.get('mkdocs', 'mkdocs_project')
    _imgs = kwargs.get('images', [])
    myMkdoc = MKDOC(os.path.join(_project_dir, _mkdocs))

    # deal with case and group
    for key, value in _project_dict.items():
        if key == 'cases':
            for _case in value:
                myMkdoc(_case, os.path.join(_project_dir, _case), append_prm=True, update=True, images=_imgs)
        else:
            my_assert(type(value) == list and value != [], TypeError, 'Input group must have a \'case\' list and it cannot be []')
            myMkdoc(key, os.path.join(_project_dir, key), append_prm=True, update=True, type='group', case_names=value, images=_imgs)
    
    # deal with analysis
    analysis_dict = kwargs.get('analysis', {})
    for key, value in analysis_dict.items():
        case_dirs = value['case_dirs']
        images = value.get('images', [])
        extra_analysis = value.get('extra_analysis', {})
        myMkdoc(key, _project_dir, append_prm=True, update=True, type='analysis', case_dirs=case_dirs, images=images, extra_analysis=extra_analysis)


def ReturnFileList(_dir, _names):
    '''
    match file starts with _name in a directory and return a list
    Inputs:
        _dir(str): target directory
        _names(list of str): names to match
    '''
    if not os.path.isdir(_dir):
        return []
    else:
        _files = []
        for _name in _names:
            _pathlist = Path(_dir).rglob(_name+'*')
            for _path in _pathlist:
                _base_name = os.path.basename(str(_path))
                if os.path.isfile(str(_path)):
                    _files.append(_base_name)
        return _files