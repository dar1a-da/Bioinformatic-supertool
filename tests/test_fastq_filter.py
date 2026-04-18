import pytest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fastq_filter import gc_percent, _normalize_bounds, mean_phred_quality, filter_fastq

def make_record(record_id: str, seq: str, phred: list[int]) -> SeqRecord:
    record = SeqRecord(Seq(seq), id=record_id, description='')
    record.letter_annotations['phred_quality'] = phred
    return record

def write_fastq(path: Path, records: list[SeqRecord]):
    with open(path, 'w') as handle:
        SeqIO.write(records, handle, 'fastq')

def test_gc_petcent_correct_value():
    assert gc_percent('GCGC') == 100.0
    assert gc_percent('ATAT') == 0.0
    assert gc_percent('ATGC') == 50.0

def test_normalize_bounds_tuple():
    assert _normalize_bounds((10, 20), default_low=0) == (10, 20)

def test_normalize_bound_value():
    assert _normalize_bounds(50, default_low=0) == (0, 50)
    assert _normalize_bounds(7.5, default_low=1.0) == (1.0, 7.5)

def test_mean_phred_quality_average():
    record = make_record('r1', 'ATGC', [10, 20, 30, 40])
    assert mean_phred_quality(record) == 25.0

def test_mean_phred_quality_zero():
    record = SeqRecord(Seq('ATGC'), id='r1', description='')
    assert mean_phred_quality(record) == 0.0

def test_fastq_filter_gc_length_quality(tmp_path):
    input_file = tmp_path / 'input.fastq'
    output_file = tmp_path / 'output.fastq'

    records = [make_record('good', 'GCGCGC', [40, 40, 40, 40, 40, 40]),
               make_record('low_gc', 'ATATAT', [40, 40, 40, 40, 40, 40]),
               make_record('short', 'GCGC', [40, 40, 40, 40]),
               make_record('low_q', 'GCGCGC', [10, 10, 10, 10, 10, 10])]
    write_fastq(input_file, records)

    result = filter_fastq(str(input_file), str(output_file),
                          gc_bounds=(40, 100), length_bounds=(5, 10),
                          quality_threshold=30)
    
    assert list(result.keys()) == ['good']


def test_fastq_filter_write_output(tmp_path):
    input_file = tmp_path / 'input.fastq'
    output_file = tmp_path / 'output.fastq'

    records = [make_record('r1', 'GCGCGC', [35, 35, 35, 35, 35, 35]),
               make_record('r2', 'ATATAT', [35, 35, 35, 35, 35, 35])]
    write_fastq(input_file, records)

    filter_fastq(str(input_file), str(output_file),
                 gc_bounds=(50, 100), length_bounds=(1, 100),
                 quality_threshold=30)
    
    assert output_file.exists()

    written_records = list(SeqIO.parse(output_file, 'fastq'))
    assert len(written_records) == 1
    assert written_records[0].id == 'r1'

def test_fastq_filter_empty(tmp_path):
    input_file = tmp_path / 'input.fastq'
    output_file = tmp_path / 'output.fastq'

    records = [make_record('r1', 'ATATAT', [10, 10, 10, 10, 10, 10]),
               make_record('r2', 'ATAT', [10, 10, 10, 10])]
    write_fastq(input_file, records)

    result = filter_fastq(str(input_file), str(output_file),
                          gc_bounds=(80, 100), length_bounds=(10, 20),
                          quality_threshold=30)
    
    assert result == {}

    written_records = list(SeqIO.parse(output_file, 'fastq'))
    assert written_records == []

def test_fastq_filter_no_file(tmp_path):
    input_file = tmp_path / 'missing.fastq'
    output_file = tmp_path / 'output.fastq'

    with pytest.raises(FileNotFoundError):
        filter_fastq(str(input_file), str(output_file),
                     gc_bounds=(0, 100), length_bounds=(0, 1000), 
                     quality_threshold=0)
        

