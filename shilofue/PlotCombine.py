# -*- coding: utf-8 -*-
r"""Combine separate figures to a bigger one.
Also combine figures from individual cases to a bigger one.

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m 

descriptions
""" 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from PIL import Image, ImageDraw, ImageFont

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    Pc_opt = PC_OPT()
    print("Combines separate figures to a bigger one.\
Also combines figures from individual cases to a bigger one.\n\
\n\
\n\
Examples of usage: \n\
\n\
  - combine figures with a json file: \n\
\n\
        Lib_PlotCombine combine_figures -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/wb_sd_issue/combine_figures.json \n\
\n\
  - options in the json file:\n\
\n\
    %s\n\
        " % Pc_opt.document())

class PC_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with json files
    List of keys:
        0: root_dir (str) - Path of the root directory
        1: case_paths (list) - relative path to individual cases
        2: plots (list) - Plots to incluce
        3: same_for_all (0 or 1) - if we use the same option of plots for all cases
        4: width (float) - width of one subplot
        5: output_dir0 (str) - Output directory
        6: if_include_case_names (0 or 1) - If we include case names
        7: title (str) - Title
        8: if_include_tile (0 or 1) - If we include a title, note this will affect all other options within the \"Title\" dict
        9: anchor (int) - The plot to anchor upon
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Path of the root directory", str, ["case_root"], '')
        self.add_key("Relative Path of case directorys", list, ["cases"], ['foo'])
        self.add_key("Plots to incluce", list, ["plots"], [])
        self.add_key("If use the same option of plots for all cases", int, ["same_for_all"], 0)
        self.add_key("Width of one subplot. This is set to -1.0 by default. By that, the width is determined with the \"anchor\" plot within the sequence",\
             float, ["width"], -1.0)
        self.add_key("Output directory", str, ["output directory"], '')
        self.add_key("If we include case names", int, ["Title", "Include case names"], 1)
        self.add_key("Title", str, ["Title", "Title"], "")
        self.add_key("If we include a title, note this will affect all other options within the \"Title\" dict",\
             int, ["Include title"], 0)
        self.add_key("The plot to anchor upon. This only takes effect when \"size\" is not given. Default is to anchor on the 0th plot in a series",\
             int, ["anchor"], 0)
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        self.n_cases = len(self.values[1])  # number of cases
        if not self.values[3]:
            # in this case, all plots need to be entered explicitly
            assert(len(self.values[2]) == self.n_cases)
            first = self.values[2][0]
            for entry in self.values[2]:
                assert(len(first) == len(entry)) # here we check they are of the same length
        pass
    
    def to_PC_init(self):
        '''
        Interface for __init__ of PC class
        '''
        root_dir = self.values[0]
        case_paths = self.values[1]
        if root_dir != '':
            for i in range(len(case_paths)):
                # substitute environmental variables and then add together
                case_paths[i] = os.path.join(Utilities.var_subs(root_dir), case_paths[i])
                Utilities.my_assert(os.path.isdir(case_paths[i]), FileNotFoundError,\
                "%s: directory %s doesn't exist" % (Utilities.func_name(), case_paths[i])) # do check
        return case_paths
    
    def to_PC_set_plots(self):
        '''
        Interface for __call__ of PC class
        '''
        if self.values[3]:
            plots = [self.values[2] for i in range(self.n_cases)]
        else:
            plots = self.values[2]
        return plots
    
    def to_PC_call(self):
        '''
        Interface for __call__ of PC class
        '''
        # size of the figure
        width = self.values[4]
        anchor = self.values[9]
        # directory to output
        root_dir = self.values[0]
        case_paths = self.values[1]
        first_case_dir = os.path.join(Utilities.var_subs(root_dir), case_paths[0])
        if self.values[5] == '':
            # in case no output dir is given, output to the first case's 'img' directory
            output_dir0 = os.path.join(first_case_dir, 'img')
        else:
            output_dir0 = Utilities.var_subs(self.values[5])
        if not os.path.isdir(output_dir0):
            os.mkdir(output_dir0)
        output_dir = os.path.join(output_dir0, 'combined')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        # title, note that None means there is no title
        if_include_title = self.values[8]
        if if_include_title:
            # title
            _title = self.values[7]
        else:
            _title = None
        # if we include a case name
        if_include_case_names = self.values[6]
        return width, anchor, output_dir, _title, if_include_case_names


class PLOT_COMBINE():
    '''
    Combine separate figures to a bigger one.
    Also combine figures from individual cases to a bigger one.
    Attributes:
        cases (list): list of case paths
        plots (2d list): plots to convert
        n_cases (int): number of cases
        n_plots (int): number of plots
        title_height (int): height of the title
    '''
    def __init__(self, case_paths):
        '''
        Initiation, read in a list of cases
        Inputs:
            case_paths (list): a list of cases
        '''
        self.cases = case_paths
        self.n_cases = len(case_paths)
        assert(self.n_cases > 0)
        self.plots = [[] for i in range(len(self.cases))]
        self.title_height = 200  # height of the title
        pass
    
    def add_case(self, case_path):
        '''
        Add one case
        Inputs:
            case_path (str): the path of a case
        '''
        self.cases.append(case_path)
        self.plots.append([])  # also add a new list of plots
        pass
    
    def set_plots(self, plots):
        '''
        Add one plot
        Inputs:
            plots (list): a list of plots
        '''
        self.plots = plots
        self.n_plots = len(plots[0])
        pass

    def configure(self):
        '''
        configuration
        '''
        pass
    
    def get_total_size(self, width, _title, **kwargs):
        '''
        get the size of new image
        Return:
            locations (list of 2 list): locations of subimages in the new combined image
            width: fixed width of a subimage
            _title (str or None) - title
        '''
        anchor = kwargs.get('anchor', 0)
        if width < 0.0:
            first_figure_path = os.path.join(self.cases[0], 'img', self.plots[0][anchor])
            Utilities.my_assert(os.path.isfile(first_figure_path), FileNotFoundError,\
             "%s: file doesn't exist (%s)" % (Utilities.func_name(), first_figure_path))
            image = Image.open(first_figure_path)
            width = image.size[0]
        locations = [[], []]
        total_size = []
        locations[0] = [i*width for i in range(self.n_cases+1)]
        if _title is None:
            locations[1].append(0)  # first one, location along the length dimension is right at the start
        else:
            locations[1].append(self.title_height)  # if there is a title, leave a little space
        for j in range(self.n_plots):
            _path = '' # initiation
            find = False
            for i in range(self.n_cases):
                # find an existing file for me
                _path = os.path.join(self.cases[i], 'img', self.plots[i][j])
                if os.path.isfile(_path):
                    find = True
                    break
            if find == False:
                # we could choose from raising an error or allow this
                # raise FileNotFoundError("No existing figure %s in all cases" % self.plots[0][j])
                length = 500  # this is the length of a vacant plot
            else:
                image = Image.open(_path)
                width_old = image.size[0]
                length_old = image.size[1]
                length = int(length_old * width / width_old) # rescale this figure
            locations[1].append(locations[1][j] + length)   # later one, add the previous location
        total_size.append(width*self.n_cases)
        total_size.append(locations[-1])
        return locations, width
        
    def draw_title(self, image, _title, if_include_case_names, w_locations):
        '''
        Draw title at the top of the image
        Inputs:
            _title (str) - title
            if_include_case_names (0 or 1) - If we include case names
            w_locations (list of int): locations along the width
        '''
        fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", 40)  # get a font
        d = ImageDraw.Draw(image)
        d.text((10,10), _title, font=fnt, fill=(0, 0, 0))  # anchor option doesn't workf
        if if_include_case_names:
            for i in range(self.n_cases):
                case_name = os.path.basename(self.cases[i])
                w_location = w_locations[i] + 10
                h_location = self.title_height / 2 + 10
                d.text((w_location,h_location), case_name, font=fnt, fill=(0, 0, 0))  # anchor option doesn't work
    
    def __call__(self, width, anchor, output_dir, _title, if_include_case_names):
        '''
        perform combination
        Inputs:
            sizes: (list of 2) - size of the plot
            output_dir: directory to output to
            _title (str or None) - title
            if_include_case_names (0 or 1) - If we include case names
        '''
        assert(os.path.isdir(output_dir))
        locations, width = self.get_total_size(width, _title, anchor=anchor)  # width is the width of a subplot
        image_size = [locations[0][-1], locations[1][-1]]
        # initiate
        new_image = Image.new('RGB',image_size,(250,250,250))
        if _title is not None:
            self.draw_title(new_image, _title,if_include_case_names, locations[0])
        for i in range(self.n_cases):
            for j in range(self.n_plots):
                plot_path = os.path.join(self.cases[i], 'img', self.plots[i][j])
                if os.path.isfile(plot_path):
                    image = Image.open(plot_path)
                    image = image.resize((width, locations[1][i+1] - locations[1][i]))  # resize to fit the spot
                else:
                    image = Image.new('RGB', (width, 500), (250,250,250)) # append a blank one
                new_image.paste(image, (locations[0][i], locations[1][j])) # paste image in place
        new_image_path = os.path.join(output_dir, 'combined_figure.png')
        print("%s: save figure: %s" % (Utilities.func_name(), new_image_path))
        new_image.save(new_image_path)
        return new_image_path


def PlotCombineFigures(json_path):
    '''
    Combine figures into a single image
    Inputs:
        json_path(str): path to the json file
    '''
    assert(os.access(json_path, os.R_OK))
    Pc_opt = PC_OPT()
    Pc_opt.read_json(json_path)  # read options
    # combine images
    Plot_Combine = PLOT_COMBINE(Pc_opt.to_PC_init())
    Plot_Combine.set_plots(Pc_opt.to_PC_set_plots())
    figure_path = Plot_Combine(*Pc_opt.to_PC_call())
    assert(os.path.isfile(figure_path))


def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='./compile_figures.json')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'combine_figures':
        # combine figures:
        assert(os.access(arg.inputs, os.R_OK))
        PlotCombineFigures(arg.inputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()