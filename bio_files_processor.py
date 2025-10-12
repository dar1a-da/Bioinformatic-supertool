import data
import os

def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=''):
    with open(input_fasta) as file:
        lines = file.readlines()
        result = []
        str = ''
        for i in range(0, len(lines)):
            if lines[i].startswith('>'):
                result.append(str+'\n')
                str = ''
                result.append(lines[i])
            else:
                str += lines[i].strip('\n')
        result.append(str)
        result.pop(0)
    
    if output_fasta != '':
        with open(os.path.join('filtered', output_fasta), mode='w') as file:
            file.writelines(result)

    return result        


print(convert_multiline_fasta_to_oneline(input_fasta='data/example_multiline_fasta.fasta'))