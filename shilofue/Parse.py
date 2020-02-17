import os 
import re
import sys

class COMPOSITION():
    """
    store value like 
      'background:4.0e-6|1.5e-6, spcrust:0.0, spharz:4.0e-6, opcrust:4.0e-6, opharz:4.0e-6 '
    or parse value back
    """
    def __init__(self, line):
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
        # Skip comment lines, mark by '#' in file
        if re.match('^ *#', line):
            pass
        # Set key and value, marked by 'set' in file
        elif re.match('^.*set', line):
            temp = re.sub('^.*set ', '', line)
            temp = temp.split('=', maxsplit=1)
            key = temp[0]
            key = key.strip(' ')
            value = temp[1]
            value = re.sub('^ *', '', value)
            value = re.sub(' *(#.*)?\n$', '', value)
            # Fix the defination of some function which
            # occupies multiple lines
            while value[-1] == '\\':
                line = fin.readline()
                line = re.sub(' *(#.*)?\n$', '', line)
                value = value + '\n' + line
            inputs[key] = value
        # Initialize new dictionary and interatively calling,
        #marked by 'subsection' in file
        elif re.match('^.*subsection', line):
            key = re.sub('^.*subsection ', '', line)
            key = key.strip('\n')
            # Fix the bug where a subsection emerges
            # multiple times
            try:
                inputs[key]
            except KeyError:
                inputs[key] = ParseFromDealiiInput(fin)
            else:
                temp = ParseFromDealiiInput(fin)
                inputs[key].update(temp.items())
        # Terminate and return, marked by 'end' in file
        elif re.match('^.*end', line):
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