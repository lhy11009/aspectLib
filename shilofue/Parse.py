import re

'''
For now, my strategy is first defining a method to parse inputs for every key word,
then defining a class to parse different type of inputs.
Todo:
    a. put them in a bigger class
    b. add method for functions
'''


class COMPOSITION():
    """
    store value like
      'background:4.0e-6|1.5e-6, spcrust:0.0, spharz:4.0e-6, opcrust:4.0e-6, opharz:4.0e-6 '
    or parse value back
    """
    def __init__(self, line):
        # parse the format:
        # key1: val1|val2, key2: val3|val4
        # to a dictionary data where
        # data[key1] = [val1, val2]
        self.data = {}
        parts = line.split(',')
        for part in parts:
            values_str = part.split(':')[1].split('|')
            # convert string to float
            values = [float(val) for val in values_str]
            self.data[part.split(':')[0]] = values

    def parse_back(self):
        """
        def parse_back(self)

        parse data back to a string
        """
        line = ''
        j = 0
        for key, values in self.data.items():
            # construct the format:
            # key1: val1|val2, key2: val3|val4
            if j > 0:
                part_of_line = ', ' + key + ':'
            else:
                part_of_line = key + ':'
            i = 0
            for val in values:
                if i == 0:
                    part_of_line += '%.4e' % val
                else:
                    part_of_line += '|' + '%.4e' % val
                i += 1
            line += part_of_line
            j += 1
        return line


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
        elif re.match('^.*set', line):
            # Parse key and value
            # from format in file as 'set key = val'
            # to a dictionary inputs
            # inputs[key] = val
            temp = re.sub('^(\t| )*set ', '', line, count=1)
            temp = temp.split('=', maxsplit=1)
            key = temp[0]
            key = key.strip(' ')
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
    return
