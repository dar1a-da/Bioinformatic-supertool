def content_gc(seq: str = 'ATGC') -> float:
    """
    Count %  GC nucleotides in sequence

    Arguments: 
    seq: str, default = 'ATGC'

    Returns: float
    """
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100


def quality_read(qual: str = 'IIII') -> float:
    """
    Calculate mean Phred quality score for a sequence

    Arguments:
    qual: str, default = 'IIII'
    
    Returns: float
    """
    phr = [ord(el) - 33 for el in qual]
    return sum(phr) / len(phr)