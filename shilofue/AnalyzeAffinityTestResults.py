r"""Analyze affinity test results

To this end, this scripts required an organized results of affinity tests to run. To run this scripts, a bash script 'organize_result.sh' should be ran beforehand.
Upon this, this scripts pull out the time results and plot figures.

This exports: 

  - a figure of results

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m shilofue.AnalyzeAffinityTestResults analyze_affinity_test_results
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_20211025 -c peloton-rome-128tasks-socket-openmpi-4.1.0
        (this includes the directory as well as the name of the cluster)

""" 
import numpy as np
import sys, os, argparse
import pathlib, subprocess
import shutil
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from shilofue.Utilities import my_assert

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']

def organize_result(test_root_dir, cluster):
    '''
    organize result
    '''
    # check directory existance
    _dir = os.path.join(test_root_dir, 'results')
    if not os.path.isdir(_dir):
        os.mkdir(_dir)
    _dir = os.path.join(test_root_dir, 'results', 'spherical_shell_expensive_solver')
    if not os.path.isdir(_dir):
        os.mkdir(_dir)
    _dir = os.path.join(test_root_dir, 'results', 'spherical_shell_expensive_solver', cluster)
    if not os.path.isdir(_dir):
        os.mkdir(_dir)
    results_dir = _dir
    results_tmp_dir = os.path.join(test_root_dir, 'tmp', cluster)
    assert(os.path.isdir(results_tmp_dir))
    for subdir, dirs, _ in os.walk(results_tmp_dir):
        for _dir in dirs:
            if _dir.startswith('output'):
                case_name = _dir.split('_', 1)[1]
                log_path = os.path.join(subdir, _dir, 'log.txt')
                target_path = os.path.join(results_dir, 'output_' + case_name)
                print("copy %s to %s" % (log_path, target_path))
                shutil.copy(log_path, target_path)
            else:
                continue

    return results_dir


def analyze_affinity_test_results(test_results_dir, output_dir):
    '''
    analyze affinity test results
    '''
    total_wall_clock = []
    assemble_stokes_system = []
    solve_stokes_system = []
    cores = []
    resolutions = []
    setups = []
    # go into sub dirs
    temp_file = os.path.join(ASPECT_LAB_DIR, 'temp')  # file to save partial results
    path_obj = pathlib.Path(test_results_dir).rglob("output*")
    i = 0
    for _path in path_obj:
        i += 1
        output_file = str(_path)
        output_path = os.path.join(test_results_dir, output_file)
        patterns = output_file.split('_')
        print("Output file found: %s" % output_path)
        # append data
        subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s" 
                       % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        try:
            data = np.genfromtxt(temp_file)
            total_wall_clock.append(data[0, -1])
            assemble_stokes_system.append(data[1, -1])
            solve_stokes_system.append(data[2, -1])
        except Exception:
            pass
        else:
            setups.append(int(patterns[-1]))
            resolutions.append(int(patterns[-2]))
            cores.append(int(patterns[-3]))
    my_assert(i > 0, AssertionError, "There is no output* file in the folder %s" % test_results_dir)

    setups = np.array(setups)
    resolutions = np.array(resolutions)
    cores = np.array(cores)
    total_wall_clock = np.array(total_wall_clock)
    assemble_stokes_system = np.array(assemble_stokes_system)
    solve_stokes_system = np.array(solve_stokes_system)
    # rearrange and sort data
    print("Affinity Test Results:")
    print('setups:', setups)
    print('resolutions', resolutions)
    print("cores:", cores)
    print("total wall clocks:", total_wall_clock)
    print("assemble_stokes_system:", assemble_stokes_system)
    print("solve_stokes_system:", solve_stokes_system)
    
    # plot via matplotlib
    resolution_options = []
    for resolution in resolutions:
        if resolution not in resolution_options:
            resolution_options.append(resolution)
    resolution_options = np.array(resolution_options)
    resolution_options = np.sort(resolution_options)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        plot_indexes = (resolutions == resolution) 
        xs = cores[plot_indexes]
        sequential_indexes = np.argsort(xs)
        ys1 = total_wall_clock[plot_indexes]
        # labels
        _label0 = 'Total Wall Clock(resolution=%d)' % resolution
        _label1 = 'Assemble Stokes System(resolution=%d)' % resolution
        _label2 = 'Solve Stokes System(resolution=%d)' % resolution
        # plot
        axs[0].loglog(xs[sequential_indexes], ys1[sequential_indexes], ".-", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label0)
        ys2 = assemble_stokes_system[plot_indexes]
        axs[0].loglog(xs[sequential_indexes], ys2[sequential_indexes], ".--", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label1)
        ys3 = ys2 / ys1
        axs[1].semilogx(xs[sequential_indexes], ys3[sequential_indexes], ".--", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label1)
        ys4 = solve_stokes_system[plot_indexes]
        axs[0].loglog(xs[sequential_indexes], ys4[sequential_indexes], ".-.", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label2)
        ys5 = ys4 / ys1
        axs[1].semilogx(xs[sequential_indexes], ys5[sequential_indexes], ".-.", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label2)
    axs[0].set_xlabel('Cores')
    axs[0].set_ylabel('Time [s]')
    axs[0].grid()
    axs[0].set_title('Wall Clock')
    axs[0].legend(fontsize='x-small')
    axs[1].set_xlabel('Cores')
    axs[1].set_ylabel('Percentage')
    axs[1].grid()
    axs[1].set_title('Percentage of Each Part')
    # title and save path 
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    plt.savefig(filepath)

    pass

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
                        default='.',
                        help='Some outputs')
    parser.add_argument('-c', '--cluster', type=str,
                        default='peloton-ii-32tasks-core-openmpi-4.0.1',
                        help='name of the cluster')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'analyze_affinity_test_results':
        # example:
        # python -m shilofue.AnalyzeAffinityTestResults analyze_affinity_test_results
        # -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_20211025 -c peloton-rome-128tasks-socket-openmpi-4.1.0
        # todo
        test_results_dir = organize_result(arg.inputs, arg.cluster)
        analyze_affinity_test_results(test_results_dir, arg.outputs)

# run script
if __name__ == '__main__':
    main()