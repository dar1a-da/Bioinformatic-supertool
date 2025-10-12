import data
from utils.module_IO_fastq import write_file
import os

def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=''):
    with open(input_fasta) as file:
        lines = file.readlines()
        result = []
        str_line = ''
        for i in range(0, len(lines)):
            if lines[i].startswith('>'):
                result.append(str_line+'\n')
                str_line = ''
                result.append(lines[i])
            else:
                str_line += lines[i].strip('\n')
        result.append(str_line)
        result.pop(0)
    
    if output_fasta != '':
        write_file(result, output_fasta)

    return result        


def parse_blast_output(input_file, output_file):
    with open(input_file) as file:
        lines = file.readlines()
        lst_description = []
        for i in range(0, len(lines)):
            if lines[i].startswith('Description'):
                point = lines[i+1].find('..')
                bracket = lines[i+1].find(']')
                if bracket == -1: bracket = point
                if point == -1: point = bracket
                index_end = min(point, bracket)
                lst_description.append(lines[i+1][:index_end+1]+'\n')

    lst_description.sort()
    write_file(lst_description, output_file)
