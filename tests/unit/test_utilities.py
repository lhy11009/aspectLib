from shilofue.Utilities import ReadHeader

def test_read_header():
    '''
    Test the ReadHeader function
    Asserts:
        key == keywords in header, with ' ' substituded by '_'
        cols == columes in header - 1
        units == units in header or None
    '''
    _texts = ['# 1: Time step', '# 2: Time (years)', '# 12: Iterations for composition solver 1']
    _header = ReadHeader(_texts)  # call function
    assert(_header['Time_step']['col'] == 0)
    assert(_header['Time_step']['unit'] == None)
    assert(_header['Time']['col'] == 1)
    assert(_header['Time']['unit'] == 'years')
    assert(_header['Iterations_for_composition_solver_1']['col'] == 11)
    assert(_header['Iterations_for_composition_solver_1']['unit'] == None)