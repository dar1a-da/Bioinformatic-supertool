from utils.module_IO_fastq import write_file

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str ='') -> list[str]:
    """
    Convert multiline fasta to oneline fasta
    
    :param input_fasta: input file fasta multiline seqs;
    :param output_fasta: output file fasta oneline seqs;
    :return: list of lines 
    """
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


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Reads a txt file, selects the first row from the Description column for each QUERY. Saves proteins to a new file in one column.

    :param input_fasta: input txt blast file;
    :param output_fasta: output txt proteins;
    """
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