import re
import os

# todo:
# a. fix this with the class that manipulates the statistic file
# b. distribute data in a separate folder

def ParseFromDealiiInput(fin):
    """
    ParseFromDealiiInput(fin)

    Parse Dealii input file to a python dictionary
    """
    inputs = {}
    line = fin.readline()
    while line is not "":
        # Inputs formats are
        # comment: "# some comment"
        # start and end of new section:
        # "subsection name" and "end"
        # set variables values:
        # 'set key = val'
        if re.match('^(\t| )*#', line):
            # Skip comment lines, mark by '#' in file
            pass
        elif re.match('^(\t| )*set', line):
            # Parse key and value
            # from format in file as 'set key = val'
            # to a dictionary inputs
            # inputs[key] = val
            temp = re.sub('^(\t| )*set ', '', line, count=1)
            temp = temp.split('=', maxsplit=1)
            key = temp[0]
            key = re.sub('(\t| )*$', '', key)
            # key = key.strip(' ')
            value = temp[1]
            value = re.sub('^ *', '', value)
            value = re.sub(' *(#.*)?\n$', '', value)
            while value[-1] == '\\':
                # Deal with entries that extent to
                # multiple lines
                line = fin.readline()
                line = re.sub(' *(#.*)?\n$', '', line)
                value = value + '\n' + line
            inputs[key] = value
        elif re.match('^.*subsection', line):
            # Start a new subsection
            # Initialize new dictionary and interatively call function,
            key = re.sub('^.*subsection ', '', line)
            key = key.strip('\n')
            try:
                # Fix the bug where a subsection emerges
                # multiple times
                inputs[key]
            except KeyError:
                inputs[key] = ParseFromDealiiInput(fin)
            else:
                temp = ParseFromDealiiInput(fin)
                inputs[key].update(temp.items())
        elif re.match('^.*end', line):
            # Terminate and return, marked by 'end' in file
            return inputs
        line = fin.readline()
    return inputs

# This class would:
# 1. Parse the variables from a prm file by saving it in a dictionary(self.idict).
# 2. Interpret from these variables to variables we give to a visit script.
# 3. Substitute the key words in a template script with these variables.
class VISIT_OPTIONS():
    """
    parse .prm file to a option file that bash can easily read
    This inherit from Utilities.CODESUB
    Attributes:
        case_dir(str): path of this case
        output_dir(str): path of the output
        visit_file(str): path of the visit file
        options(dict): dictionary of key and value to output
    """
    def __init__(self, case_dir):
        """
        Initiation
        Args:
            case_dir(str): directory of case
        """
        # check directory
        self.case_dir = case_dir
        assert(os.path.isdir(self.case_dir))
        self.output_dir = os.path.join(case_dir, 'output')
        assert(os.path.isdir(self.output_dir))
        self.visit_file = os.path.join(self.output_dir, 'solution.visit')
        assert(os.access(self.visit_file, os.R_OK))
        # output dir
        self.output_dir = os.path.join(case_dir, 'output')
        if not os.path.isdir(self.output_dir):  # mkdir
            os.mkdir(self.output_dir)
        # img dir
        self.img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(self.img_dir):
            os.mkdir(self.img_dir)
        # get inputs from .prm file, assume prm file is in the case directory
        prm_file = os.path.join(self.case_dir, 'case.prm')
        assert(os.access(prm_file, os.R_OK))
        with open(prm_file, 'r') as fin:
            self.idict = ParseFromDealiiInput(fin)
        # initiate a dictionary
        self.options = {}
        # initiate a statistic data
        # self.Statistics = Plot.LINEARPLOT('Statistics')
        # self.statistic_file = os.path.join(self._output_dir, 'statistics')
        # self.Statistics.ReadHeader(self.statistic_file)
        # self.Statistics.ReadData(self.statistic_file)

        # horiz_avg
        # self.horiz_avg_file = os.path.join(self._output_dir, "depth_average.txt")
    
    def Interpret(self, **kwargs):
        """
        Interpret the inputs parsed from a prm file
        kwargs: options
            last_steps(list): plot the last few steps
        """
        # visit file
        self.options["VISIT_FILE"] = self.visit_file
        # particle file
        particle_file = os.path.join(self.output_dir, 'particles.visit')
        if os.access(particle_file, os.R_OK):
            self.options["VISIT_PARTICLE_FILE"] = particle_file
        # directory to output data
        self.options["DATA_OUTPUT_DIR"] = self.output_dir
        # directory to output images
        if not os.path.isdir(self.img_dir):
            os.mkdir(self.img_dir)
        self.options["IMG_OUTPUT_DIR"] = self.img_dir
        # own implementations
        # initial adaptive refinement
        self.options['INITIAL_ADAPTIVE_REFINEMENT'] = self.idict['Mesh refinement'].get('Initial adaptive refinement', '6')
        # get snaps for plots
        graphical_snaps, _, _ = GetSnapsSteps(self.case_dir, 'graphical')
        self.options['ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS'] = str(graphical_snaps)
        # plot slab 
        self.options['IF_PLOT_SLAB'] = 'True'
        last_step = graphical_snaps[-1] - int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])
        last_steps = kwargs.get('last_steps', None)
        if type(last_steps) == int:
            self.options['PLOT_SLAB_STEPS'] = [i for i in range(last_step - last_steps + 1, last_step + 1)]
        else:
            self.options['PLOT_SLAB_STEPS'] = [0, 1, 2, 3, 4, 5, 6, 7]
        # self.options['IF_DEFORM_MECHANISM'] = value.get('deform_mechanism', 0)
        self.options['IF_DEFORM_MECHANISM'] = 1

# todo: fix this with the class that manipulates the statistic file
def GetSnapsSteps(case_dir, type_='graphical'):
    case_output_dir = os.path.join(case_dir, 'output')
    # parse parameters from case.prm
    prm_file = os.path.join(case_dir, 'case.prm')
    assert(os.access(prm_file, os.R_OK))
    with open(prm_file, 'r') as fin:
        idict = ParseFromDealiiInput(fin)
    # import statistics file
    Statistics = Plot.STATISTICS_PLOT_OLD('Statistics')
    statistic_file = os.path.join(case_output_dir, 'statistics')
    Utilities.my_assert(os.access(statistic_file, os.R_OK), FileNotFoundError,
              'case statistic file - %s cannot be read' % prm_file)
    Statistics.ReadHeader(statistic_file)
    Statistics.ReadData(statistic_file)
    col_time = Statistics.header['Time']['col']

    # final time
    final_time = Statistics.data[-1, col_time]

    # time interval
    # graphical
    try:
        time_between_graphical_output = float(idict['Postprocess']['Visualization']['Time between graphical output'])
    except KeyError:
        time_between_graphical_output = 1e8
    total_graphical_outputs = int(final_time / time_between_graphical_output) + 1
    graphical_times = [i*time_between_graphical_output for i in range(total_graphical_outputs)]
    graphical_steps = [Statistics.GetStep(time) for time in graphical_times]
    # particle
    try:
        time_between_particles_output = float(idict['Postprocess']['Particles']['Time between data output'])
    except KeyError:
        time_between_particles_output = 1e8
    total_particles_outputs = int(final_time / time_between_particles_output) + 1
    particle_times = [i*time_between_particles_output for i in range(total_particles_outputs)]
    particle_steps = [Statistics.GetStep(time) for time in particle_times]
    # initial_snap
    try:
        initial_snap = int(idict['Mesh refinement']['Initial adaptive refinement'])
    except KeyError:
        initial_snap = 6
    # end snap
    snaps = [0]
    if type_ == 'graphical':
        start_ = initial_snap
        end_ = total_graphical_outputs + initial_snap
        snaps = list(range(start_, end_))
        times = graphical_times
        steps = graphical_steps
    elif type_ == 'particle':
        start_ = 0
        end_ = total_particles_outputs
        snaps = list(range(start_, end_))
        times = particle_times
        steps = particle_steps
    return snaps, times, steps

# Initiate the class, Parse the variables from a prm file(case.prm) by saving it in a dictionary(self.idict).
visit_options = VISIT_OPTIONS('.')
# visit_options.Interpret(last_steps=3)