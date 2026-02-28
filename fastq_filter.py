from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord


def gc_percent(seq: str) -> float:
    return gc_fraction(seq) * 100

def _normalize_bounds(bounds, *, default_low):
    if isinstance(bounds, tuple):
        low, high = bounds
    else:
        low, high = default_low, bounds
    return low, high

def mean_phred_quality(record: SeqRecord) -> float:
    """Mean quality of read by Phred"""
    quals = record.letter_annotations.get('phred_quality')
    if not quals:
        return 0.0
    return sum(quals) / len(quals)


def filter_fastq(input_fastq: str, output_fastq:str,
                 gc_bounds: tuple[float, float] | float = (0,100),
                 length_bounds: tuple[int, int] | float = (0,2**32),
                 quality_threshold: float = 0) -> dict[str, SeqRecord]:
    
    """FASTQ filtering by GC (%), length, average Phred quality.
    Returns dict: {record.id : SeqRecord} for those who passed the filter.
    Writes the past reads to output_fastq (fastq format)."""

    gc_low, gc_high = _normalize_bounds(gc_bounds, default_low=0.0)
    len_low, len_high = _normalize_bounds(length_bounds, default_low=0)

    passed: list[SeqRecord] = []
    result = dict[str, SeqRecord] = {}

    for record in SeqIO.parse(input_fastq, 'fastq'):
        seq_str = str(record.seq)
        L = len(seq_str)
        gc = gc_percent(seq_str)
        q = mean_phred_quality(record)

        gc_ok = gc_low <= gc <= gc_high
        len_ok = len_low <= L <= len_high
        q_ok = q >= quality_threshold

        if gc_ok and len_ok and q_ok:
            passed.append(record)
            result[record.id] = record

    with open(output_fastq, 'w') as out_handle:
        SeqIO.write(passed, out_handle, 'fastq')

    return result
