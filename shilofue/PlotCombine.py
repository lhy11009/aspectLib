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
import json
import shutil
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import patches as mpatches
from PIL import Image, ImageDraw, ImageFont
from shilofue.PlotRunTime import PlotFigure as RunTimePlotFigure
from shilofue.PlotStatistics import STATISTICS_PLOT

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("Combines separate figures to a bigger one.\
Also combines figures from individual cases to a bigger one.\n\
\n\
\n\
Examples of usage: \n\
\n\
  - combine figures with a json file (from different cases): \n\
\n\
        Lib_PlotCombine combine_figures -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/wb_sd_issue/combine_figures.json \n\
\n\
  - prepare_results using the IMAGE_OPT module (from one case, e.g. combine plots from one step)\n\
        Lib_PlotCombine prepare_results -i ~/ASPECT_PROJECT/aspectLib/files/TwoDSubduction/211230/figure_option.json \n\
\n\
  - prepare results by combining runtime outputs\n\
        Lib_PlotCombine combine_runtime -j `pwd`/runtime.json \n\
\n\
        ")


def ShowJsonOption():
    Pc_opt = PC_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Pc_opt.document_str()
        )

####
# A base class for the task of combining plots
####
class PC_OPT_BASE(Utilities.JSON_OPT):
    def __init__(self):
        '''
        initiation of the class
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Path of the root directory", str, ["case_root"], '', nick="case_root")
        self.add_key("Relative Path of case directorys", list, ["cases"], ['foo'], nick="case_relative_paths")
        self.add_key("Width of one subplot. This is set to -1.0 by default. By that, the width is determined with the \"anchor\" plot within the sequence",\
             float, ["width"], -1.0, nick="width")
        self.add_key("Use relative path for the output directory", int, ["output directory", "relative"],\
        0, nick="output_dir_relative")
        self.add_key("Path for the output directory", str, ["output directory", "path"], '', nick="output_dir_path")
        self.n_cases = len(self.values[1])  # number of cases
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        # check existence of cases
        case_root = self.values[0]
        case_relative_paths = self.values[1]
        for case_relative_path in case_relative_paths:
            case_path = os.path.join(Utilities.var_subs(case_root), case_relative_path)
            if not os.path.isdir(case_path):
                raise FileExistsError("Directory %s doesn't exist." % case_path)
        # make the output directory if not existing
        output_dir_relative = self.values[3]
        output_dir_path = self.values[4]
        if output_dir_relative == 1:
            output_dir = os.path.join(Utilities.var_subs(case_root), output_dir_path)
        else:
            output_dir = Utilities.var_subs(output_dir_path)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    
    def get_case_absolute_paths(self):
        '''
        return the absolute paths of cases
        '''
        case_root = self.values[0]
        case_relative_paths = self.values[1]
        case_absolute_paths = []
        for case_relative_path in case_relative_paths:
            case_path = os.path.join(Utilities.var_subs(case_root), case_relative_path)
            case_absolute_paths.append(case_path)
        return case_absolute_paths

    def get_case_names(self):
        '''
        return the names of cases
        '''
        case_names = []
        for case_relative_path in case_relative_paths:
            case_name = os.path.basename(case_relative_path)
            case_names.append(case_name)
        return case_names

    def get_output_dir(self):
        '''
        return the output dir
        '''
        case_root = self.values[0]
        case_relative_paths = self.values[1]
        output_dir_relative = self.values[3]
        output_dir_path = self.values[4]
        if output_dir_relative == 1:
            output_dir = os.path.join(Utilities.var_subs(case_root), output_dir_path)
        else:
            output_dir = Utilities.var_subs(output_dir_path)
        return output_dir
    
    def get_color_json_output_path(self):
        '''
        return a path to output the color scheme to a json file
        '''
        output_dir = self.get_output_dir()
        color_json_output_path = os.path.join(output_dir, 'color.json')
        return color_json_output_path


####
# Classes and functions for the previous function - combine figures from different figures in a big one
####

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
        self.add_key("Name of the plot", str, ["Name"], "foo", nick="_name")
    
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
        # name of the plot
        _name = self.values[10]
        return width, anchor, output_dir, _title, if_include_case_names, _name


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
                This is where the upper-left corner of each figure is located. Index 0 in
                the first dimension is along the horizontal direction and index 1 is along the
                vertical direction.
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
        print("locations: ", locations) # debug
        return locations, width
        
    def draw_title(self, image, _title, if_include_case_names, w_locations):
        '''
        Draw title at the top of the image
        Inputs:
            _title (str) - title
            if_include_case_names (0 or 1) - If we include case names
            w_locations (list of int): locations along the width
        '''
        individual_figure_width = w_locations[1] - w_locations[0]
        if False:
            fnt_size = 40
        else:
            fnt_size = 40
        fnt0 = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", fnt_size)  # get a font
        d = ImageDraw.Draw(image)
        d.text((10,10), _title, font=fnt0, fill=(0, 0, 0))  # anchor option doesn't workf
        if if_include_case_names:
            for i in range(self.n_cases):
                case_name = os.path.basename(self.cases[i])
                str_len = fnt0.getsize(case_name)
                if str_len[0] > 0.9 * individual_figure_width:
                    # resize font with longer case name
                    fnt_size_re = int(fnt_size * 0.9 * individual_figure_width // str_len[0])
                    fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", fnt_size_re)
                else:
                    fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", fnt_size)
                w_location = w_locations[i] + 10
                h_location = self.title_height / 2 + 10
                d.text((w_location,h_location), case_name, font=fnt, fill=(0, 0, 0))  # anchor option doesn't work
    
    def __call__(self, width, anchor, output_dir, _title, if_include_case_names, _name):
        '''
        perform combination
        Inputs:
            sizes: (list of 2) - size of the plot
            output_dir: directory to output to
            _title (str or None) - title
            if_include_case_names (0 or 1) - If we include case names
            _name (str) - name of the plot
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
                    image = image.resize((width, locations[1][j+1] - locations[1][j]))  # resize to fit the spot
                else:
                    image = Image.new('RGB', (width, 500), (250,250,250)) # append a blank one
                new_image.paste(image, (locations[0][i], locations[1][j])) # paste image in place
        new_image_path = os.path.join(output_dir, '%s.png' % _name)
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


def PrepareResults(json_path):
    '''
    Prepare results using the IMAGE_OPT module.
    Inputs:
        json_path(str): path to the json file
    '''
    assert(os.path.isfile(json_path))
    pillow_opt = Utilities.PILLOW_OPT()
    pillow_opt.read_json(json_path)
    Utilities.PillowRun(*pillow_opt.to_pillow_run())

####
# 06292022: new functions combines one kind of figures at a time (e.g. run time, visualization)
####
class PC_RUNTIME_OPT(PC_OPT_BASE):
    '''
    Define a class to work with json files
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        PC_OPT_BASE.__init__(self)
        self.start = self.number_of_keys()

    def to_init(self):
        '''
        interfaces to the __init__ function
        '''
        case_absolute_paths = self.get_case_absolute_paths()
        return case_absolute_paths

    def to_call(self):
        '''
        interfaces to the __call__ function
        '''
        width = self.values[2]
        output_dir = self.get_output_dir()
        return width, output_dir



class PLOT_COMBINE_RUNTIME(PLOT_COMBINE):
    '''
    Combine results from run time outputs
    '''
    def __init__(self, case_paths):
        PLOT_COMBINE.__init__(self, case_paths)
        UnitConvert = Utilities.UNITCONVERT()
        self.StatisticPlotter = STATISTICS_PLOT("statistic_plot", unit_convert=UnitConvert)
        pass

    def __call__(self, width, output_dir, **kwargs):
        '''
        perform combination
        Inputs:
            sizes: (list of 2) - size of the plot
            output_dir: directory to output to
        '''
        _name = "combine_runtime"
        _title = "Comparing run-time results"
        n_color_max = 10
        color_method = kwargs.get('color_method', 'generated')
        dump_color_to_json = kwargs.get('dump_color_to_json', None)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        # initiate
        case_names = []  # names of cases
        for i in range(self.n_cases):
            case_name = os.path.basename(self.cases[i])
            case_names.append(case_name)
        ni = 4  # number of plots along 1st and 2nd dimension
        nj = 2
        fig = plt.figure(tight_layout=True, figsize=[5*nj, 5*ni])
        gs = gridspec.GridSpec(ni, nj)
        colors_dict = {}
        if color_method == 'generated':
            colors_dict['max'] = n_color_max
            if self.n_cases > n_color_max:
                raise ValueError("max number of colors must be bigger than the number of cases")
            normalizer = [ float(i)/(n_color_max) for i in range(self.n_cases) ]
            colors = cm.rainbow(normalizer)
            for i in range(self.n_cases):
                case_name = os.path.basename(self.cases[i])
                colors_dict[case_name] = list(colors[i])
        elif color_method == 'check_first':
            colors = []
            if dump_color_to_json is not None:
                if os.path.isfile(dump_color_to_json):
                    with open(dump_color_to_json, 'r') as fin:
                        colors_dict = json.load(fin)
                else:
                    colors_dict['max'] = n_color_max
            # first loop to get the number of colors in the json file
            n_color_in_json = 0
            for i in range(self.n_cases):
                case_name = os.path.basename(self.cases[i])
                if case_name in colors_dict:
                    n_color_in_json += 1
            normalizer = [ float(i)/(n_color_max) for i in range(n_color_in_json, self.n_cases) ]
            new_colors = cm.rainbow(normalizer)
            # second loop to assign colors to new cases
            j = 0
            for i in range(self.n_cases):
                case_name = os.path.basename(self.cases[i])
                try:
                    colors.append(colors_dict[case_name])
                except KeyError:
                    colors.append(new_colors[j])
                    colors_dict[case_name] = list(colors[i])
                    j += 1
        else:
            raise ValueError('Not implemented')
        # plot number of cells
        ax = fig.add_subplot(gs[1, 0])
        for i in range(self.n_cases):
            case_name = os.path.basename(self.cases[i])
            # plot results and combine
            statistic_file_path = os.path.join(self.cases[i], 'output', 'statistics')
            self.StatisticPlotter.ReadData(statistic_file_path)
            self.StatisticPlotter.ReadHeader(statistic_file_path)
            self.StatisticPlotter.PlotNumberOfCells(axis=ax, color=colors[i])
            pass
        ax.legend()
        # plot degree of freedoms
        ax = fig.add_subplot(gs[1, 1])
        for i in range(self.n_cases):
            if i == 0:
                label_all = True
            else:
                label_all = False
            case_name = os.path.basename(self.cases[i])
            # plot results and combine
            statistic_file_path = os.path.join(self.cases[i], 'output', 'statistics')
            self.StatisticPlotter.ReadData(statistic_file_path)
            self.StatisticPlotter.ReadHeader(statistic_file_path)
            self.StatisticPlotter.PlotDegreeOfFreedom(axis=ax, color=colors[i], label_all=label_all)
            pass
        ax.legend()
        # plot temperatures
        ax = fig.add_subplot(gs[2, 0])
        for i in range(self.n_cases):
            if i == 0:
                label_all = True
            else:
                label_all = False
            case_name = os.path.basename(self.cases[i])
            # plot results and combine
            statistic_file_path = os.path.join(self.cases[i], 'output', 'statistics')
            self.StatisticPlotter.ReadData(statistic_file_path)
            self.StatisticPlotter.ReadHeader(statistic_file_path)
            self.StatisticPlotter.PlotTemperature(axis=ax, color=colors[i], label_all=label_all)
            pass
        ax.legend()
        # plot velocities
        ax = fig.add_subplot(gs[2, 1])
        for i in range(self.n_cases):
            if i == 0:
                label_all = True
            else:
                label_all = False
            case_name = os.path.basename(self.cases[i])
            # plot results and combine
            statistic_file_path = os.path.join(self.cases[i], 'output', 'statistics')
            self.StatisticPlotter.ReadData(statistic_file_path)
            self.StatisticPlotter.ReadHeader(statistic_file_path)
            self.StatisticPlotter.PlotVelocity(axis=ax, color=colors[i], label_all=label_all)
            pass
        ax.legend()
        # plot number of non-linear iterations
        ax = fig.add_subplot(gs[3, 0])
        for i in range(self.n_cases):
            if i == 0:
                label_all = True
            else:
                label_all = False
            case_name = os.path.basename(self.cases[i])
            # plot results and combine
            statistic_file_path = os.path.join(self.cases[i], 'output', 'statistics')
            self.StatisticPlotter.ReadData(statistic_file_path)
            self.StatisticPlotter.ReadHeader(statistic_file_path)
            self.StatisticPlotter.PlotNumberOfNonlinearIterations(axis=ax, color=colors[i], label_all=label_all)
            pass
        ax.legend()
        # plot run time info
        ax = fig.add_subplot(gs[0, 1])
        ax_twin = ax.twinx()
        for i in range(self.n_cases):
            if i == 0:
                label_all = True
                append_extra_label = True
                if_legend = True
            else:
                label_all = False
                append_extra_label = False
                if_legend = False
            case_name = os.path.basename(self.cases[i])
            # plot results and combine
            log_file_path = os.path.join(self.cases[i], 'output', 'log.txt')
            RunTimePlotFigure(log_file_path, None, savefig=False, axis=ax, fix_restart=True,\
            label_all=label_all, append_extra_label=append_extra_label, if_legend=if_legend,\
            twin_axis=ax_twin, color=colors[i], x_variable='time')
            pass
        # plot the color labels
        ax = fig.add_subplot(gs[0, 0])
        PlotColorLabels(ax, case_names, colors)
        # generate figures
        fig_path = os.path.join(output_dir, '%s.png' % _name)
        print("%s: save figure: %s" % (Utilities.func_name(), fig_path))
        plt.savefig(fig_path)
        # dump a color file
        if dump_color_to_json is not None:
            assert(os.path.isdir(os.path.dirname(dump_color_to_json)))
            with open(dump_color_to_json, 'w') as fout:
                json.dump(colors_dict, fout)
            print("%s: dump color options: %s" % (Utilities.func_name(), dump_color_to_json))
        return fig_path

def PlotColorLabels(ax, case_names, colors): 
    '''
    plot the color labels used for different cases
    '''
    labels = []
    patches = []
    for case_name in case_names:
        labels.append(case_name)
    for _color in colors:
        patches.append(mpatches.Patch(color=_color))
    ax.legend(patches, labels, loc='center', frameon=False)
    ax.axis("off")


def PlotCombineRuntime(json_path):
    '''
    Combine runtime plot
    Inputs:
        json_path: path to a json file for configuration
    '''
    assert(os.access(json_path, os.R_OK))
    Pc_opt = PC_RUNTIME_OPT()
    Pc_opt.read_json(json_path)  # read options
    # plot the combined figure
    PlotCombineRunTime = PLOT_COMBINE_RUNTIME(Pc_opt.to_init())
    PlotCombineRunTime(*Pc_opt.to_call(), dump_color_to_json=Pc_opt.get_color_json_output_path(),\
    color_method='check_first')
    # save the configuration file
    json_copy_path = os.path.join(Pc_opt.get_output_dir(), 'runtime.json')
    try:
        shutil.copy(json_path, json_copy_path)
        print("Saved json file: ", json_copy_path)
    except shutil.SameFileError:
        print("Saing json file: The two files are the same, skip")
    


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
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='A json file for configuration')
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
    elif _commend == 'prepare_results':
        # prepare results
        assert(os.access(arg.inputs, os.R_OK))
        PrepareResults(arg.inputs)
    elif _commend == "combine_runtime":
        PlotCombineRuntime(arg.json)
    elif (_commend in ['--json_option', '-jo']):
        # json options
        ShowJsonOption()
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()