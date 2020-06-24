from shilofue.Doc import SeparateNavPattern, ExtractNav, ProduceNav
from shilofue.Utilities import re_count_indent


def test_seperate_nav_pattern():
    '''
    test the function SeparateNavPattern(_pattern)
    '''
    # test 1
    _pattern = ' - test: foo  \n'
    key, value = SeparateNavPattern(_pattern)
    assert(key == 'test' and value == 'foo')
    # test 2
    _pattern = ' \t- test:  \t\n'
    key, value = SeparateNavPattern(_pattern)
    assert(key == 'test' and value == '')


def test_extract_nav():
    '''
    test the function ExtractNav(_lines, **kwargs)
    '''
    # test 1
    _lines = [' - Home: index.md',
              ' - About:',
              '   - About: about/about.md',
              '   - About1: about/about1.md',
              ' - Foo: foo.md']
    _odict, _at= ExtractNav(_lines, previous=re_count_indent(_lines[0]))
    assert(_at == 4)
    assert(_odict == {'Home': 'index.md', 
                      'About': {'About': 'about/about.md',
                                'About1': 'about/about1.md'},
                      'Foo': 'foo.md'})
    # test 2
    _lines = [' - Home: index.md',
              ' - About:',
              '   - About: about/about.md',
              '   - About1: about/about1.md',
              '   - About2:',
              '        - About2: about/about2.md',
              '        - About3: about/about3.md',
              ' - Foo: foo.md']
    _odict, _at= ExtractNav(_lines, previous=re_count_indent(_lines[0]))
    assert(_at == 7)
    assert(_odict == {'Home': 'index.md', 
                      'About': {'About': 'about/about.md',
                                'About1': 'about/about1.md',
                                'About2': {'About2': 'about/about2.md',
                                           'About3': 'about/about3.md'}},
                      'Foo': 'foo.md'})


def test_produce_nav():
    '''
    test the function ProduceNav(_lines, **kwargs)
    '''
    _idict = {'Home': 'index.md', 
              'About': {'About': 'about/about.md',
              'About1': 'about/about1.md',
              'About2': {'About2': 'about/about2.md',
                         'About3': 'about/about3.md'}},
              'Foo': 'foo.md'}
    _lines_standard = \
             ['    - Home: index.md',
              '    - About:',
              '        - About: about/about.md',
              '        - About1: about/about1.md',
              '        - About2:',
              '            - About2: about/about2.md',
              '            - About3: about/about3.md',
              '    - Foo: foo.md']
    _lines = ProduceNav(_idict)
    print('lines: ', _lines)  # screen output
    assert(_lines==_lines_standard)