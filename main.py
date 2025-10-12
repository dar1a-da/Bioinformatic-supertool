from data import seqs
from utils.module_filter_fastq import content_gc, quality_read, compare
from utils.module_dna_rna_tools import check_nucleic_acid, transcribe_nucleic_acid, reverse_nucleic_acid, complement_nucleic_acid, reverse_complement_nucleic_acid

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
    for idname, (seq, quality) in seqs.items():
        gc_res = compare(gc_bounds,seq,content_gc)
        len_res = compare(length_bounds,seq,len)

        if quality_threshold is not None:
            phr_res = quality_threshold <= quality_read(quality)
        else:
            phr_res == True

        if gc_res and len_res and phr_res:
            result[idname] = seq, quality
    
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
    
    operations = {
    "is_nucleic_acid": check_nucleic_acid,
    "transcribe": transcribe_nucleic_acid,
    "reverse": reverse_nucleic_acid,
    "complement": complement_nucleic_acid,
    "reverse_complement": reverse_complement_nucleic_acid
    }
    
    results = []
    for seq in sequences:
        func = operations.get(process)
        if process == "is_nucleic_acid":
            results.append(func(seq))
        elif func and check_nucleic_acid(seq):
            results.append(func(seq))
        
    if len(results) == 1:
        return results[0]
    return results