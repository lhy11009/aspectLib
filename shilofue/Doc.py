import numpy
import json
import os
import re
import shilofue.json
from importlib import resources
from shutil import copyfile
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
        case_name(str) - case name
        dir(str) - directory of case data
        imgs(list) - list of images
        new_files(dict) - name: new file
    '''
    def __init__(self, _case_name, _dir):
        '''
        Inputs:
            _case_name(str) - case name
            _dir(str) - directory of case data
        '''
        self.case_name = _case_name
        self.dir = _dir
        self.imgs = []
        self.new_files = {}

    def append_img(self, _img):
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

    def __call__(self, _odir, **kwargs):
        '''
        Call function and write to file that will be used bu mkdoc
        Inputs:
            kwargs:
                update(bool) - if update an existing case
                append_prm(bool) - if append a prm file
                prm(str) - name of a prm file
        '''
        update = kwargs.get('update', False)
        append_prm = kwargs.get('append_prm', False)
        _prm = kwargs.get('prm', 'case.prm')
        # These directory and files need to be preexist
        assert(os.path.isdir(_odir))
        assert(os.path.isdir(os.path.join(_odir, 'docs')))
        _mcdocs_file = os.path.join(_odir, 'mkdocs.yml')
        assert(os.path.isfile(_mcdocs_file))
        _index_file = os.path.join(_odir, 'docs', 'index.md')
        _case_dir = os.path.join(_odir, 'docs', self.case_name)
        if not os.path.isdir(_case_dir):
            os.mkdir(_case_dir)
        # Append a summary.md
        if os.path.isfile(os.path.join(_case_dir, 'summary.md')):
            if update == True:
                self.GenerateCaseMkd(_case_dir, update=True)
        else:
            _filename = self.GenerateCaseMkd(_case_dir)
            # in a mkdocs file, files are listed as 'case_name/_filename'
            self.new_files['Summary'] = os.path.join(self.case_name, os.path.basename(_filename))
        # Append a prm file if the append_prm option is True
        if append_prm:
            if os.path.isfile(os.path.join(_case_dir, _prm)):
                if update == True:
                    _prm_source_file = os.path.join(self.dir, _prm)
                    assert(os.path.isfile(_prm_source_file))
                    copyfile(_prm_source_file, os.path.join(_case_dir, _prm))
            else:
                _prm_source_file = os.path.join(self.dir, _prm)
                assert(os.path.isfile(_prm_source_file))
                copyfile(_prm_source_file, os.path.join(_case_dir, _prm))
                self.new_files['Parameters'] = os.path.join(self.case_name, _prm)
        if self.new_files != {}:
            self.RenewMkdocsYml(_odir)

    def GenerateCaseMkd(self, _case_dir, **kwargs):
        '''
        Generate markdown file of a case
        Inputs:
            _case_dir(str): directory of this case
            kwargs:
                filename(str): name of the file
        Returns:
            _filename(str): file generated
        '''
        _filename = kwargs.get('filename', 'summary.md')
        _filename = os.path.join(_case_dir, _filename)
        _auto_mkd_file = os.path.join(self.dir, 'auto.md')
        _extra_mkd_file = os.path.join(self.dir, 'extra.md')
        assert (os.path.isfile(_auto_mkd_file))
        copyfile(_auto_mkd_file, _filename)  # copy the auto.mkd file from case directory
        if (os.path.isfile(_extra_mkd_file)):
            # apppend contents of extra.mkd at the end of case.mkd
            with open(_filename, 'a') as fout:
                with open(_extra_mkd_file, 'r') as fin:
                    fout.write(fin.read())
        return _filename

    def RenewMkdocsYml(self, _odir):
        '''
        Renew the mkdocs.yml file in the project directory
        Inputs:
            _odir(str): project directory
        '''
        _filename = os.path.join(_odir, 'mkdocs.yml')
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
            value = _nav_dict[self.case_name]
            my_assert(type(value) == dict, TypeError,
                      'entry for a case must be a single dict, prepared to include dictionary in the future')
            value = {**value, **self.new_files}  # merge and substitute value in first dict with value in second dict
        except KeyError:
            # this case is new to the yml file
            _nav_dict[self.case_name] = self.new_files
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
                _odict[key], i = ExtractNav(_lines, previous=re_count_indent(_lines[i+1]), at=i+1)
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
    for key, value in _project_dict.items():
        if key == 'cases':
            for _case in value:
                myMkdoc = MKDOC(_case, os.path.join(_project_dir, _case))
                myMkdoc(os.path.join(_project_dir, _mkdocs), append_prm=True, update=True)
