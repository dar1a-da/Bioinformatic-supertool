from data import seqs
from utils.module_filter_fastq import *
from utils.module_dna_rna_tools import *

def filter_fastq(
        seqs: dict[str, tuple[str, str]] = {
            'id1': ('ATGC', 'IIII')
        }, 
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
    result = {}
    for key, (seq, quality) in seqs.items():
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
            result[key] = seq, quality
    
    return result



def run_dna_rna_tools(*args):
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
