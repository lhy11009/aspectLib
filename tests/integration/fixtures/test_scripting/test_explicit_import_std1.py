def PlotFigure(file_path, fig_path):
    '''
    descriptions
    Inputs:
        - file_path(str): path of a statistic file of aspect
        - figure_path(str): path of the output figure
    Returns:
        -
    '''
    # Init the UnitConvert class
    UnitConvert = Utilities.UNITCONVERT()
    json_file = os.path.join(ASPECT_LAB_DIR, 'shilofue', 'json_files', 'post_process.json')
    with open(json_file, 'r') as fin:
        pdict = json.load(fin)
    plot_options = pdict.get('Statistics', {})
    Plotter = STATISTICS_PLOT('Statistics', unit_convert=UnitConvert, options=plot_options)
    fig_dir = os.path.dirname(fig_path)
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
    fig_generated_path = Plotter(file_path, fileout=fig_path)  # plot figure
    print("New figure: %s" % fig_generated_path)
    return fig_generated_path 
    pass
