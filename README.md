# Bioinformatic-supertool

This repository contains tools for working with nucleotide sequences.

**The entry point** to the program is a file `main.py `, which implements calls to the main functions.

The script `bio_files_processor.py` in which the functions `convert_multiline_fasta_to_oneline`, `parse_blast_output` are realized.

- In folder `data` file with fasta sequences (`example_multiline_fasta.fasta`), a file with fastq format sequences: name, sequences, "+" string, string quality (`example_fastq.fastq`), file with blast result (`example_blast_result.txt`), file (`example_gbk.gbk`).
- In folder `utils` there are modules with additional functions (`module_*.py`).
- In folder `filtered` function results are stored `filter_fastq` (**filt**), `convert_multiline_fasta_to_oneline` (**one_line.fasta**), `parse_blast_output` (**proteins.txt**).

## Function `dna_rna_tools`

Function allows to perform basic operations with the DNA/RNA sequence:
- check whether the string is a nucleotide sequence;
- return the transcribed, reverse, complementary, reverse complementary sequence.

The **dna_rna_tools** function accepts an arbitrary number of nucleotide sequences (`str`) and the type of process to be performed on the sequences (`str`),
`process` — type of operation (`"is_nucleic_acid"`, `"transcribe"`, `"reverse"`, `"complement"`, `"reverse_complement"`).

**Usage Example**

```python
run_dna_rna_tools('ATGC', 'is_nucleic_acid') # True
run_dna_rna_tools('AGU', 'CCuu' 'transcribe') # ['AGT', 'CCtt']
run_dna_rna_tools('ATg', 'reverse') # 'gTA'
run_dna_rna_tools('ctA', 'complement') # 'gaT'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
```

## Function `filter_fastq`

Function works with sequences in the **fastq** format. Allows to filter them by GC composition, reading length, and read quality on the Phred33 scale.

The **filter_fastq** function accepts 5 arguments as input:
1. The name of the file with fastq format sequences.
2. The name of the file where the result will be recorded.
3. Interval (tuple of two values) or GC value (`float') of the composition (in percent) for filtering. In the case of a single border, all reads below this border are filtered (saved).
4. The interval (tuple of two values) or the value (`float`) of the reed length for filtering. In the case of a single border, all reads below that border are filtered.
5. The threshold value of the average read quality for filtering (`int`). Reads with quality below the threshold value are discarded.

Saves fastq sequence strings filtered by the set parameters.

```python
filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0)
```

**Usage Example**
```python
filter_fastq('input_fastq', 'output_fastq', gc_bounds = (0,100),        length_bounds = (0,2**32), quality_threshold = 0)
```

## Function `convert_multiline_fasta_to_oneline`

Function reads the input fasta file, in which the sequence (DNA/RNA/protein) can be split into several lines, and then saves it to a new fasta file in which each sequence fits into one line.

Accepts 2 arguments as input (input_fasta and output_fasta). output_fasta optional argument, if not specified, returns a list.

```python
convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str ='') -> list[str]
```

## Function `parse_blast_output`

The program reads the preset txt file, selects the first row from the Description column for each QUERY (paragraph Sequences producing significant alignments:). The set of obtained proteins is saved in alphabetical order to a new file.

Accepts 2 arguments as input (input_file, output_file). 

```python
parse_blast_output(input_file: str, output_file: str) -> None
```
