import re
import os

# In the following function, the inputs(prm file) is parsed to a dictionary
# Each subsection in the prm file would be parsed to sub-dictionary relative to its upper level
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


# Compare the contents of 'inputs' below with the prm file to see what it is.
# the 'Dimension' variable locates at the top level and it is a single entry in the dictionary.
# 'Solver parameters' is a name for a subsection and it is a sub-dictionary itself.
prm_path = 'case.prm'
assert(os.access(prm_path, os.R_OK))
with open(prm_path, 'r') as fin:
    inputs = ParseFromDealiiInput(fin)
print('Contents of inputs:\n', 'Dimension:', inputs['Dimension'], '\n', 'Solver parameters:\n', inputs['Solver parameters'])


# In the following function, a dictionary is paresed back to inputs(prm file)
def ParseToDealiiInput(fout, outputs, layer=0):
    """
    def ParseToDealiiInput(fout, outputs, layer=0)

    Parse a python dictionary into a Dealii input file
    """
    indent = ' ' * 4 * layer  # Indentation of output
    for key, value in outputs.items():
        # output format is
        # same as input but no comment
        if type(value) is str:
            fout.write(indent + 'set %s = %s\n' % (key, value))
        elif type(value) is dict:
            if layer == 0:
                fout.write('\n')
            fout.write(indent + 'subsection %s\n' % key)
            layer1 = layer + 1
            ParseToDealiiInput(fout, value, layer1)
            fout.write(indent + 'end\n')
            if layer == 0:
                fout.write('\n')
        else:
            raise ValueError('Value in dict must be str')
    return


# Next we would change an option (i.e. lower the level of refinement by 1) in the prm file and output it to another prm file
# Check the case_o.prm file to see if this works
outputs = inputs.copy()
outputs['Mesh refinement']['Initial global refinement'] = '4'
prm_path_o = 'case_o.prm'
with open(prm_path_o, 'w') as fout:
    ParseToDealiiInput(fout, outputs)


