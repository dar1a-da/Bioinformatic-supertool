import argparse
import logging
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord


logging.basicConfig(filename='fastq_filter.log', level=logging.INFO, filemode='a')

logger = logging.getLogger(__name__)

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

    logger.info('Start filtering')

    passed: list[SeqRecord] = []
    result: dict[str, SeqRecord] = {}

    for record in SeqIO.parse(input_fastq, 'fastq'):
        seq_str = str(record.seq)
        length = len(seq_str)
        gc = gc_percent(seq_str)
        q = mean_phred_quality(record)

        gc_ok = gc_low <= gc <= gc_high
        len_ok = len_low <= length <= len_high
        q_ok = q >= quality_threshold

        if gc_ok and len_ok and q_ok:
            passed.append(record)
            result[record.id] = record

    with open(output_fastq, 'w') as out_handle:
        SeqIO.write(passed, out_handle, 'fastq')

    logger.info('Success')

    return result


def parse_args():
    parser = argparse.ArgumentParser(description='FASTQ filtering by GC (%), length, average Phred quality')

    parser.add_argument('input_fastq', help='Path to input Fastq file')

    parser.add_argument('output_fastq', help='Path to output Fastq file')

    parser.add_argument('--gc-bounds', nargs=2, type=float, metavar=('GC_min', 'GC_max'),
                        default=(0.0, 100.0), help='Bounds of procent GC, f.e. --gc-bounds 40 60')
    
    parser.add_argument('--length-bounds', nargs=2, type=int, metavar=('LEN_min', 'LEN_max'), 
                        default=(0, 2**32), help='Length bounds, f.e. --length-bounds 50 150')

    parser.add_argument('--quality-threshold', type=float, default=0.0, 
                        help='Average quality of read')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    try:
        result = filter_fastq(input_fastq=args.input_fastq, output_fastq=args.output_fastq,
                          gc_bounds=tuple(args.gc_bounds), length_bounds=tuple(args.length_bounds),
                          quality_threshold=args.quality_threshold)
    
        print(f"Passed reads: {len(result)}")
    except Exception as e:
        logger.error('Failed', exc_info=True)
        print(f"Error {e}")


#python3 fastq_filter.py data/example_fastq.fastq filtered/output.fastq --gc-bounds 40 60 --length-bounds 50 150 --quality-threshold 30
#Passed reads: 22
