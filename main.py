#from data import seqs
import os
from utils.module_filter_fastq import *
from utils.module_dna_rna_tools import *
from utils.module_IO_fastq import *

def filter_fastq(
        input_fastq,
        output_fastq,
        gc_bounds: tuple[float, float] | float = (0,100),
        length_bounds: tuple[int, int] | float = (0,2**32),
        quality_threshold: float = 0
        ) -> dict[str, tuple[str, str]]:
    """
    Filter reads by gc composition, length, Phred quality scale

    Arguments:
    seqs: dict: {'id': ('ATGC', 'IIII'}, default={'id1': ('ATGC', 'IIII')}
    gc_bounds: tuple or float: (min, max) or number, default=(0, 100)
    length_bounds: tuple or int: (min, max) or number, default=(0,2**32)
    quality_threshold: float, default=0

    Returns:
    dict[str, tuple[str, str]]
    dict filter sequences
    """
    seqs = read_fastq(input_fastq)
    result = {}
    for key, (seq, plus, quality) in seqs.items():
        if isinstance(gc_bounds, tuple):
            gc_res = (gc_bounds[0] <= content_gc(seq)) and (gc_bounds[1] >= content_gc(seq))
        else:
            gc_res = content_gc(seq) <= gc_bounds

        if isinstance(length_bounds, tuple):
            len_res = (length_bounds[0] <= len(seq)) and (length_bounds[1] >= len(seq))
        else:
            len_res = len(seq) <= length_bounds

        if quality_threshold is not None:
            phr_res = quality_threshold <= quality_read(quality)
        else:
            phr_res == True

        if gc_res and len_res and phr_res:
            result[key] = seq, plus, quality
    
    write_fastq(result, output_fastq)
    return result


def run_dna_rna_tools(*args: str) -> bool | list[bool] | str | list[str]:
    """
    Perform basic DNA/RNA sequence operations.

    This function check if strings are nucleotide sequences, 
    return their transcribed, reverse, complementary, or reverse-complementary forms.

    Parameters:
    *args : str
        Variable number of arguments.
        All arguments except the last one are nucleotide sequences (str),
        the last argument specifies the operation:
        - 'is_nucleic_acid'
        - 'transcribe'
        - 'reverse'
        - 'complement'
        - 'reverse_complement'

    Returns:
    bool or list[bool] or str or list[str]
        - If process is 'is_nucleic_acid': returns bool or list[bool].
        - Otherwise: returns the resulting sequence(s) as str or list[str].
    """
    *sequences, process = args
    results = []
    
    for seq in sequences:
        if process == 'is_nucleic_acid':
            results.append(check_nucleic_acid(seq))

        elif process == 'transcribe' and check_nucleic_acid(seq):
            results.append(transcribe_nucleic_acid(seq))

        elif process == 'reverse' and check_nucleic_acid(seq):
            results.append(reverse_nucleic_acid(seq))

        elif process == 'complement' and check_nucleic_acid(seq):
            results.append(complement_nucleic_acid(seq))

        elif process == 'reverse_complement' and check_nucleic_acid(seq):
            results.append(reverse_complement_nucleic_acid(seq))

    if len(results) == 1:
        return results[0]
    return results



filter_fastq(input_fastq='data/example_fastq2.fastq', output_fastq='filt')

