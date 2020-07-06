import numpy
import json
import os
import re
import shilofue.json
from importlib import resources
from shutil import copyfile
from pathlib import Path
from shilofue.Utilities import my_assert, re_neat_word, re_count_indent



class DOC():
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
            with resources.open_text(shilofue.json, _layout_json) as fin:
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
        self.imgs = kwargs.get('images', [])
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
        '''
        self.new_files = {}
        update = kwargs.get('update', False)
        append_prm = kwargs.get('append_prm', False)
        _prm = kwargs.get('prm', 'case.prm')
        _type = kwargs.get('type', 'case')
        _case_names = kwargs.get('case_names', None)
        # These directory and files need to be preexist
        assert(os.path.isdir(self.odir))
        assert(os.path.isdir(os.path.join(self.odir, 'docs')))
        _mcdocs_file = os.path.join(self.odir, 'mkdocs.yml')
        assert(os.path.isfile(_mcdocs_file))
        _index_file = os.path.join(self.odir, 'docs', 'index.md')
        if _type == 'case':
            _target_dir = os.path.join(self.odir, 'docs', _name)
            self.AppendCase(_name, _dir, _target_dir, update=update, append_prm=append_prm)
        elif _type == 'group':
            my_assert(_case_names is not None, ValueError, 'For a group, case names cannot be None. Valid names must be given')
            _target_dir = os.path.join(self.odir, 'docs', _name)
            self.AppendGroup(_name, _dir, _case_names, _target_dir, update=update, append_prm=append_prm)
        else:
            raise ValueError('Type must be \'case\' or \'group\'')
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
            self.new_files[_case_name] = {}
            _case_dir = os.path.join(_dir, _case_name)
            _case_target_dir = os.path.join(_target_dir, _case_name)
            self.AppendCase(_case_name, _case_dir, _case_target_dir, update=update, append_prm=append_prm, basename=_name)

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
        assert (os.path.isfile(_auto_mkd_file))
        copyfile(_auto_mkd_file, _filename)  # copy the auto.mkd file from case directory
        if (os.path.isfile(_extra_mkd_file)):
            # apppend contents of extra.mkd at the end of case.mkd
            with open(_extra_mkd_file, 'r') as fin:
                _contents = fin.read()
            with open(_filename, 'a') as fout:
                fout.write('\n')
                fout.write(_contents)
        if _imgs != []:
            # append images
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
            value = {**value, **self.new_files}  # merge and substitute value in first dict with value in second dict
        except KeyError:
            # this case is new to the yml file
            _nav_dict[_name] = self.new_files
        _new_lines = _lines[0: _start]
        _new_lines += ProduceNav(_nav_dict)
        _new_lines += _lines[_end: -1]
        with open(_filename, 'w') as fout:
            for _line in _new_lines:
                fout.write(_line + '\n')


def ExtractNav(_lines, **kwargs):
    '''
    extract the nav message from plain text
    Inputs:
        _contents(str): plain text contents
        kwargs:
            todo
        todo:
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
    _temp = _pattern.split('-')[1]
    key = re_neat_word(_temp.split(':')[0])
    value = re_neat_word(_temp.split(':')[1])
    return key, value



def UpdateProjectDoc(_project_dict, _project_dir, **kwargs):
    '''
    Update doc for all cases in this project
    '''
    _mkdocs = kwargs.get('mkdocs', 'mkdocs_project')
    _imgs = kwargs.get('images', [])
    myMkdoc = MKDOC(os.path.join(_project_dir, _mkdocs), images=_imgs)
    for key, value in _project_dict.items():
        if key == 'cases':
            for _case in value:
                myMkdoc(_case, os.path.join(_project_dir, _case), append_prm=True, update=True)
        else:
            my_assert(type(value) == list and value != [], TypeError, 'Input group must have a \'case\' list and it cannot be []')
            myMkdoc(key, os.path.join(_project_dir, key), append_prm=True, update=True, type='group', case_names=value)


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
                _files.append(_base_name)
        return _files