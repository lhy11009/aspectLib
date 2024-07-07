# -*- coding: utf-8 -*-
r"""Plot run time information

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage:

        python -m shilofue.PlotRunTime plot
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/output/log.txt
        -o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/img/run_time.png

  - plot newton solver history
         python -m shilofue.PlotRunTime plot_newton_solver_history
         -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20
         -o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20/img/newton_solver_history.png

descriptions
"""
from re import S
import numpy as np
import sys, os, argparse
import subprocess
import numpy as np
import shilofue.Plot as Plot
import warnings
from matplotlib import pyplot as plt
from matplotlib import gridspec
import pdb

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
This scripts generate plots\n\
\n\
Examples of usage: \n\
\n\
  - plot newton solver results\n\
    -s: assign a minimum step; \n\
    -s1: assign a maximum step \n\
\n\
        Lib_PlotRunTime plot_newton_solver_history -i `pwd` -s 350 -s1 450\n\
\n\
  - plot newton solver step\n\
\n\
        Lib_PlotRunTime plot_newton_solver_step -i `pwd` -s 390\n\
        ")

def RunTimeInfo(log_path, **kwargs):
    '''
    Read runtime info from log file and then print them
    Inputs:
        log_path(str) - path to log file
    '''
    save_temp_file_local = kwargs.get('save_temp_file_local', False)
    quiet = kwargs.get("quiet")
    if save_temp_file_local:
        temp_dir = os.path.join(os.path.dirname(log_path), 'py_outputs')
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
        temp_path = os.path.join(temp_dir, 'run_time_output')
    else:
        temp_path = os.path.join(RESULT_DIR, 'run_time_output')
    trailer = None # add this to filename
    hr = 3600.0  # hr to s
    # read log file
    if not quiet:
        print("awk -f %s/bash_scripts/awk_states/parse_block_output %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_output %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))

    Plotter = Plot.LINEARPLOT('RunTime', {})
    Plotter.ReadHeader(temp_path)
    Plotter.ReadData(temp_path)

    # give warning if there is no data
    # future: reorder error messages
    if Plotter.data.size == 0:
        warning_message = "func %s: There is no block data in file %s, Do nothing" % (Utilities.func_name(), log_path)
        warnings.warn(warning_message, Utilities.WarningTypes.FileHasNoContentWarning)
        return 0

    # fix indexes
    col_time = Plotter.header['Time']['col']
    unit_time = Plotter.header['Time']['unit']
    col_step = Plotter.header['Time_step_number']['col']
    col_wallclock = Plotter.header['Wall_Clock']['col']
    unit_wallclock = Plotter.header['Wall_Clock']['unit']
    times = Plotter.data[:, col_time]
    steps = Plotter.data[:, col_step]
    wallclocks = Plotter.data[:, col_wallclock]
    
    # last step info
    # todo
    last_step = steps[-1]
    last_time = times[-1] 
    last_wallclock = wallclocks[-1]
    if not quiet:
        print("step: %d, time: %.4e, wallclock: %.4e" % (last_step, last_time, last_wallclock))
    return last_step, last_time, last_wallclock

def PlotFigure(log_path, fig_path, **kwargs):
    '''
    Read runtime info from log file and then plot
    Inputs:
        log_path(str) - path to log file
    Returns:
        RunTime.png
    '''
    trailer = None # add this to filename
    to_myr = 1e6
    fix_restart = kwargs.get('fix_restart', False)
    savefig = kwargs.get('savefig', True)
    _color = kwargs.get('color', None)
    ax1 = kwargs.get('axis', None)
    ax1_twinx = kwargs.get('twin_axis', None)
    ax2 = kwargs.get('ax2', None) # wallclock - step gradient
    ax3 = kwargs.get('ax3', None) # wallclock - time gradient
    if_legend = kwargs.get('if_legend', False)
    x_variable = kwargs.get('x_variable', 'step')
    save_temp_file_local = kwargs.get('save_temp_file_local', False)  # save temp file to a local place
    if save_temp_file_local:
        temp_dir = os.path.join(os.path.dirname(log_path), 'py_outputs')
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
        temp_path = os.path.join(temp_dir, 'run_time_output')
    else:
        temp_path = os.path.join(RESULT_DIR, 'run_time_output')
    if ax1 == None:
        fig, ax1 = plt.subplots(figsize=(5, 5))
    append_extra_label = kwargs.get('append_extra_label', True)  # add this because multiple plots mess up the axis
    hr = 3600.0  # hr to s
    # read log file
    print("awk -f %s/bash_scripts/awk_states/parse_block_output %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_output %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))

    Plotter = Plot.LINEARPLOT('RunTime', {})
    Plotter.ReadHeader(temp_path)
    Plotter.ReadData(temp_path)

    # give warning if there is no data
    # future: reorder error messages
    if Plotter.data.size == 0:
        warning_message = "func %s: There is no block data in file %s, Do nothing" % (Utilities.func_name(), log_path)
        warnings.warn(warning_message, Utilities.WarningTypes.FileHasNoContentWarning)
        return 0

    # fix indexes
    col_time = Plotter.header['Time']['col']
    unit_time = Plotter.header['Time']['unit']
    col_step = Plotter.header['Time_step_number']['col']
    col_wallclock = Plotter.header['Wall_Clock']['col']
    unit_wallclock = Plotter.header['Wall_Clock']['unit']
    times = Plotter.data[:, col_time]
    steps = Plotter.data[:, col_step]
    wallclocks = Plotter.data[:, col_wallclock]


    # fix bugs in data
    i = 1
    n_found = 0
    n_size = times.shape[0]
    while i < n_size - n_found:
        wallclock_pre = wallclocks[i-1]
        step_pre = steps[i-1]
        wallclock = wallclocks[i]
        step = steps[i]
        if step > step_pre and wallclock < wallclock_pre :
            # this is an error, delete this one
            times[i: n_size-1] = times[i+1: n_size]
            steps[i: n_size-1] = steps[i+1: n_size]
            wallclocks[i: n_size-1] = wallclocks[i+1: n_size]
            n_found += 1
        else:
            i += 1
    times = times[0: n_size-n_found]
    steps = steps[0: n_size-n_found]
    wallclocks = wallclocks[0: n_size-n_found]

    # fix restart
    re_inds = []
    steps_fixed = np.array([])  # initialize these 3 as the original ones, so we'll see no changes if there is no restart
    times_fixed = np.array([])
    wallclocks_fixed = np.array([])
    if fix_restart:
        last_step = -1
        i = 0
        for step in steps:
            if step <= last_step:
                re_inds.append(i)  # step < last step is a marker for restart
            last_step = step
            i = i+1
        if re_inds != []:
            for i in range(len(re_inds)-1):
                re_ind = re_inds[i]
                re_ind_next = re_inds[i+1]
                wallclocks[re_ind: re_ind_next] += wallclocks[re_ind - 1]
            re_ind = re_inds[-1]  # deal with the last one seperately
            wallclocks[re_ind: ] += wallclocks[re_ind - 1]
            # steps_fixed = steps[re_inds]
            steps_fixed = steps
            # times_fixed = times[re_inds]
            times_fixed = times
            # wallclocks_fixed = wallclocks[re_inds]
            wallclocks_fixed = wallclocks
        else:
            steps_fixed = steps
            times_fixed = times
            wallclocks_fixed = wallclocks

    # mask for time
    t_mask = (times >= 0.0)  # should always be true
    t_mask_fixed = (times_fixed >= 0.0)
    try:
        time_range = kwargs['time_range']
    except KeyError:
        pass
    else:
        Utilities.my_assert(((type(time_range) == list) and (len(time_range) == 2)), TypeError,\
                  "PlotFigure: time_range should be a list of 2")
        t_mask = ((times >= time_range[0]) & (times <= time_range[1]))
        t_mask_fixed = ((times_fixed >= time_range[0]) & (times_fixed <= time_range[1]))
        trailer = "%.2e_%.2e" % (time_range[0], time_range[1])

    # line 1: time
    # use mask
    line_type = '-'
    if _color is None:
        color = 'tab:blue'
        color_tab = 'tab:blue'
        line_type = '-'
    else:
        color = _color
        color_tab = 'k'
        line_type = '-'
    if x_variable == 'step':
        ax1.plot(steps[t_mask], times[t_mask] / to_myr, line_type, color=color, label='Time')
        ax1.plot(steps_fixed[t_mask_fixed], times_fixed[t_mask_fixed] / to_myr, '.', color=color)
        ax1.set_xlabel('Step')
        ax1.set_ylabel('Time [myr]', color=color_tab)
    elif x_variable == 'time':
        ax1.plot(times[t_mask] / to_myr, steps[t_mask], line_type, color=color, label='Step')
        ax1.plot(times_fixed[t_mask_fixed] / to_myr, steps_fixed[t_mask_fixed], '.', color=color)
        ax1.set_xlabel('Time [myr]', color=color_tab)
        ax1.set_ylabel('Step')
    else:
        raise ValueError("x_variable needs to be \"time\" or \"step\"")
    ax1.tick_params(axis='y', labelcolor=color_tab)
    ax1.set_title('Run Time')
    if if_legend:
        ax1.legend()
    # line 2: wall clock
    if ax1_twinx == None:
        ax1_twinx = ax1.twinx()
    if _color is None:
        color = 'tab:red'
        color_tab = 'tab:red'
        line_type = '-'
    else:
        color = _color
        color_tab = 'k'
        line_type = '--'
    if x_variable == 'step':
        ax1_twinx.plot(steps[t_mask], wallclocks[t_mask] / hr, line_type, color=color, label='Wall Clock [s]')
        ax1_twinx.plot(steps_fixed[t_mask_fixed], wallclocks_fixed[t_mask_fixed] / hr, '.', color=color)
    elif x_variable == 'time':
        ax1_twinx.plot(times[t_mask] / to_myr, wallclocks[t_mask] / hr, line_type, color=color, label='Wall Clock [s]')
        ax1_twinx.plot(times_fixed[t_mask_fixed] / to_myr, wallclocks_fixed[t_mask_fixed] / hr, '.', color=color)
    else:
        raise ValueError("x_variable needs to be \"time\" or \"step\"")  
    if append_extra_label:
        ax1_twinx.set_ylabel('Wall Clock [hr]', color=color_tab)
        if x_variable == 'step':
            ax1_twinx.set_xlabel('Step')
        elif x_variable == 'time':
            ax1_twinx.set_xlabel('Time [myr]')
        ax1_twinx.tick_params(axis='y', labelcolor=color_tab)
    if if_legend:
        ax1_twinx.legend()
    
    # plot wallclock / step gradient
    if ax2 is not None:
        assert(fix_restart is True) # only work if the restart steps are fixed
        # additional plot: time gradient
        ax2_twinx = ax2.twinx()
        wc_gradient_fixed = np.gradient(wallclocks_fixed/hr, steps_fixed)
        t_gradient_fixed = np.gradient(times_fixed / to_myr, steps_fixed)
        ax2.plot(steps_fixed, t_gradient_fixed, '--', color='tab:blue')
        ax2_twinx.plot(steps_fixed, wc_gradient_fixed, '--', color='tab:red')
        ax2.set_xlabel('Step')
        ax2.set_ylabel('Time Gradient [myr / step]', color='tab:blue')
        ax2_twinx.set_ylabel("Wall Clock Gradient [hr / step]", color='tab:red')
        ax2.tick_params(axis='y', labelcolor="tab:blue")
        ax2_twinx.tick_params(axis='y', labelcolor="tab:red")

    # plot wallclock / time gradient
    if ax3 is not None:
        assert(fix_restart is True) # only work if the restart steps are fixed
        # additional plot: time gradient
        gradient_fixed = np.gradient(wallclocks_fixed/hr, times_fixed/to_myr)
        ax3.plot(times_fixed / to_myr, gradient_fixed, '--', color='tab:red')
        ax3.grid()
        ax3.set_xlabel('Time [myr]')
        ax3.set_ylabel("Wall Clock Gradient [hr / myr]", color='tab:red')
        ax3.tick_params(axis='y', labelcolor="tab:red")

    # save figure
    if savefig: 
        _name, _extension = Utilities.get_name_and_extention(fig_path)
        if trailer != None:
            fig_path = "%s_%s.%s" % (_name, trailer, _extension)
        else:
            fig_path = fig_path
        fig.tight_layout()
        plt.savefig(fig_path)
        print("New figure: %s" % fig_path)
    else:
        fig_path=None
    # return range of steps
    steps_plotted = steps[t_mask]
    return fig_path, [steps_plotted[0], steps_plotted[-1]]


def PlotNewtonSolver(log_path, fig_path_base, **kwargs):
    '''
    Read runtime info from log file and then plot
    Inputs:
        log_path(str) - path to log file
    Returns:
        RunTime.png
    '''
    # read log file
    temp_path = os.path.join(RESULT_DIR, 'run_time_output_newton')
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_newton %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))
    Plotter = Plot.LINEARPLOT('RunTime', {})
    Plotter.ReadHeader(temp_path)
    Plotter.ReadData(temp_path)

    step = kwargs.get('step', 0)  # step to plot
    col_step = Plotter.header['Time_step_number']['col']
    mask_step = (Plotter.data[:, col_step] == step)
    data = Plotter.data[mask_step, :]

    col_number_of_iteration = Plotter.header['Index_of_nonlinear_iteration']['col']
    number_of_iterations = data[:, col_number_of_iteration]
    col_residual = Plotter.header['Relative_nonlinear_residual']['col']
    residuals = data[:, col_residual]
    col_scaling_factor = Plotter.header["Newton_Derivative_Scaling_Factor"]['col']
    scaling_factors = data[:, col_scaling_factor]
    # line1: residual
    fig, ax1 = plt.subplots(figsize=(5, 5))
    color = 'tab:blue'
    ax1.semilogy(number_of_iterations, residuals, '-', color=color, label='Residuals')
    ax1.set_ylabel('Relative non-linear residual', color=color)
    ax1.set_xlabel('Number of non-linear iterations')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_title('Newton Solver Output For Step %d' % (step))
    # line 2: scaling factor
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.plot(number_of_iterations, scaling_factors, '--', color=color, label='Scaling factors')
    ax2.set_ylabel('Newton derivative scaling factor', color=color)
    ax2.set_xlabel('Number of non-linear iterations')
    ax2.tick_params(axis='y', labelcolor=color)

    # save figure
    fig.tight_layout()
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    fig_path = "%s_s%06d.%s" % (fig_path_base0, step, fig_path_type)
    plt.savefig(fig_path)
    print("New figure: %s" % fig_path)
    return fig_path


def PlotNewtonSolverHistory(log_path, fig_path_base, **kwargs):
    '''
    Read runtime info from log file and then plot
    Inputs:
        log_path(str) - path to log file
    Returns:
        SolverHistory.png
    '''
    # read log file
    temp_path = os.path.join(RESULT_DIR, 'run_time_output_newton')
    
    query_iterations = kwargs.get('query_iterations', None)  # number of iteration to query
    
    trailer = None
    
    print("awk -f %s/bash_scripts/awk_states/parse_block_newton %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_newton %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))
    Plotter = Plot.LINEARPLOT('SolverHistory', {})
    Plotter.ReadHeader(temp_path)
    try:
        Plotter.ReadData(temp_path)
    except ValueError as e:
        raise ValueError('Value error(columns are not uniform): check file %s' % temp_path)
    col_step = Plotter.header['Time_step_number']['col']
    col_number_of_iteration = Plotter.header['Index_of_nonlinear_iteration']['col']
    col_residual = Plotter.header['Relative_nonlinear_residual']['col']
    end_step = int(Plotter.data[-1, col_step])
    
    # steps to plot, todo get from default or read in
    steps = np.array([i for i in range(end_step)])
    number_of_iterations = np.zeros(end_step)
    residuals = np.zeros(end_step)

    # outputs for querying additional iterations 
    n_query_iteration = 0
    if query_iterations is not None:
        for iteration in query_iterations:
            assert(type(iteration) == int)
        n_query_iteration = len(query_iterations)
        residuals_at_iterations = np.zeros([n_query_iteration, end_step])

    for i in range(steps.size):
        step = steps[i]
        mask_step = (Plotter.data[:, col_step] == step)
        data = Plotter.data[mask_step, :]
        number_of_iterations[step] = data[-1, col_number_of_iteration]
        residuals[step] = data[-1, col_residual]
        
    # query residual of iteration 'query_iteration'
    for i in range(steps.size):
        step = steps[i]
        mask_step = (Plotter.data[:, col_step] == step)
        data = Plotter.data[mask_step, :]
        for j in range(n_query_iteration):
            try:
                residuals_at_iterations[j, step] = data[query_iterations[j], col_residual]
            except IndexError:
                residuals_at_iterations[j, step] = data[-1, col_residual]

    # plot mask
    step_range = kwargs.get('step_range', None)
    if step_range == None:
        s_mask = (steps >= 0)
    else:
        Utilities.my_assert(type(step_range) == list and len(step_range) == 2, TypeError, "%s: step_range must be a list of 2." % Utilities.func_name())
        s_mask = ((steps >= step_range[0]) & (steps <= step_range[1]))  # this is hard coded to be 0 for now
        trailer = "%d_%d" % (step_range[0], step_range[1])

    # todo_newton
    # line1: residual
    fig = plt.figure(tight_layout=True, figsize=(15, 5))  # plot of wallclock
    gs = gridspec.GridSpec(1, 3)
    color = 'tab:blue'

    # figure 1: plot the linear results, residual and scaling factor
    # residual of iteration, plot mask
    ax = fig.add_subplot(gs[0, 0])
    ax.semilogy(steps[s_mask], residuals[s_mask], '-', linewidth=1.5, color=color, label='Residuals')
    # query residual of iteration 'query_iteration'
    if query_iterations is not None:
        for j in range(n_query_iteration):
            ax.semilogy(steps[s_mask], residuals_at_iterations[j, s_mask], '--', linewidth=0.5, label='Residuals at iteration %d' % query_iterations[j])
    ax.set_ylabel('Relative non-linear residual', color=color)
    ax.set_xlabel('Steps')
    ax.set_title('Solver History')
    ax.tick_params(axis='y', labelcolor=color)
    ax.legend()
    # line 2: scaling factor
    ax2 = ax.twinx()
    color = 'tab:red'
    ax2.plot(steps[s_mask], number_of_iterations[s_mask], '.', color=color, label='Numbers of Iterations')
    ax2.set_ylabel('Numbers of Iterations', color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    # figure 2: plot the linear results, residual and scaling factor
    # bin plot of the logrithm value of residuals.
    # bin size is twiked to 0.25
    ax = fig.add_subplot(gs[0, 1])
    data = np.log10(residuals[s_mask])
    plt.hist(data, bins=np.arange(int(np.floor(min(data))), int(np.floor(max(data))) + 1, 0.25), edgecolor='black')
    # Add labels and title
    plt.xlabel('log(relative non-linear residuals)')
    plt.ylabel('Frequency')
    plt.title('Bin Plot')

    # figure 2: plot the linear results, residual and scaling factor
    # bin plot of the residulas
    # bin size is twiked to 10
    ax = fig.add_subplot(gs[0, 2])
    data = number_of_iterations[s_mask]
    plt.hist(data, bins=np.arange(0, int(np.ceil(max(data)/10.0 + 1.0)*10.0), 10.0), color='red', edgecolor='black')
    # Add labels and title
    plt.xlabel('Numbers of Iterations')
    plt.ylabel('Frequency')
    plt.title('Bin Plot')

    # save figure
    fig.tight_layout()
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    if trailer == None:
        fig_path = "%s.%s" % (fig_path_base0, fig_path_type)
    else:
        fig_path = "%s_%s.%s" % (fig_path_base0, trailer, fig_path_type)
    plt.savefig(fig_path)
    print("New figure (new): %s" % fig_path)
    return fig_path


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
                        help='Some inputs')
    parser.add_argument('-o', '--outputs', type=str,
                        default='',
                        help='Some outputs')
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='time step')
    parser.add_argument('-s1', '--step1', type=int,
                        default=-1,
                        help='time step')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot':
        # example:
        PlotFigure(arg.inputs, arg.outputs, fix_restart=True)

    elif _commend == 'plot_case':
        # example:
        log_file = os.path.join(arg.inputs, 'output', 'log.txt')
        assert(log_file)
        img_dir = os.path.join(arg.inputs, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        fig_path = os.path.join(img_dir, 'run_time.png')

        fig = plt.figure(tight_layout=True, figsize=(6, 18))
        gs = gridspec.GridSpec(3,1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        PlotFigure(log_file, fig_path, axis=ax1, ax2=ax2, ax3=ax3, savefig=False, fix_restart=True) 
        fig.savefig(fig_path)
    
        print("New figure: %s" % fig_path)
    
    elif _commend == 'case_run_time_info':
        # example:
        log_file = os.path.join(arg.inputs, 'output', 'log.txt')
        assert(log_file)
        # todo
        RunTimeInfo(log_file)

    elif _commend == "plot_newton_solver_step":
        # plot newton solver output
        log_file = os.path.join(arg.inputs, 'output', 'log.txt')
        assert(log_file)
        fig_path = os.path.join(arg.inputs, 'img', 'newton_solver.png')
        if not os.path.isdir(os.path.dirname(fig_path)):
            os.mkdir(os.path.dirname(fig_path))
        o_path = PlotNewtonSolver(log_file, fig_path, step=arg.step)

    elif _commend == "plot_newton_solver_history":
        # plot newton solver output
        log_file = os.path.join(arg.inputs, 'output', 'log.txt')
        assert(log_file)
        fig_path = os.path.join(arg.inputs, 'img', 'newton_solver_history.png')
        if not os.path.isdir(os.path.dirname(fig_path)):
            os.mkdir(os.path.dirname(fig_path))
        if arg.step1 >= 0:
            o_path = PlotNewtonSolverHistory(log_file, fig_path, step_range=[arg.step, arg.step1])
        else:
            o_path = PlotNewtonSolverHistory(log_file, fig_path)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()
    else:
        # commend not right
        raise ValueError('%s is not a valid commend' % _commend)

# run script
if __name__ == '__main__':
    main()
