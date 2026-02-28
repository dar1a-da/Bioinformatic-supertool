# Bioinformatic-supertool

This repository contains tools for working with nucleotide sequences.

## Content
- Requirements
- Classes
  - InvalidAlphabetError
  - BiologicalSequence
  - NucleicAcidSequence
  - DNASequence
  - RNASequence
  - AminoAcidSequence
- Function filter_fastq
- Usage examples

## Requirements
- Python 3.10+
- Biopython 1.86

## Classes
### InvalidAlphabetError
A custom exception (the successor of ValueError) that is thrown when invalid characters are found in the sequence.
Occurs:
When creating a sequence object (DNASequence("ATBX")),
When you call check_alphabet()

### BiologicalSequence
An abstract base class that captures a common interface.
Records:
len(obj) — the length of the sequence  
indexing obj[i] — character by index  
slices of obj[i:j] — returns an object of the same class  
beautiful output via __str__ and __repr__  
check_alphabet() method (required for heirs)  
Fields:  
seq: str — string with sequence  
ALPHABET: ClassVar[frozenset[str]] — valid alphabet (set in the heirs)

### NucleicAcidSequence
The base class for nucleic acids. The heir of BiologicalSequence.  
Implements:
check_alphabet() — alphabet check (uses ALPHABET heir),  
complement() — complementary sequence (polymorphically via COMPLEMENT_MAP),  
reverse() — the reverse sequence,  
reverse_complement() — reverse complement.  
Important:  
Polymorphism is implemented via the COMPLEMENT_MAP class attribute.  
The complement() method does not contain if DNA/RNA conditions, and returns an object of the desired type.

### DNASequence
The successor of NucleicAcidSequence for DNA.
transcribe() — transcribes DNA into RNA and returns RNASequence.

### RNASequence
The successor of NucleicAcidSequence for RNA.

### AminoAcidSequence
A class for amino acid sequences. The heir of BiologicalSequence.
Implements:  
check_alphabet() — checking the amino acid alphabet  
translate_aa(nuc, frame=0, stop_at_stop=True, unknown="X") is a class-method for translating DNA/RNA into protein  
Broadcast Parameters:  
frame: 0/1/2 — offset of the reading frame  
stop_at_stop: if True, stop at the first stop code. *  
unknown: character for an unknown codon (default is "X")  
Returns: the AminoAcidSequence object.

## Function filter_fastq (Biopython)
Filters reads in FASTQ by:  
length (length_bounds)  
average Phred quality (quality_threshold)  
GC-composition (gc_bounds)  
Implemented via:  
Bio.SeqIO.parse, Bio.SeqIO.write  
Bio.SeqUtils.gc_fraction  
qualities are taken from record.letter_annotations["phred_quality"]

## Usage examples
```python
dna = DNASequence('ATGCCGTA')
print(dna)                 # DNASequence(len=8, seq='ATGCCGTA')
print(dna[0])              # 'A'
print(dna[2:6])            # DNASequence(len=4, seq='GCCG')

dna = DNASequence('ATGC')
print(dna.complement())    # DNASequence(len=4, seq='TACG')
print(dna.reverse())       # DNASequence(len=4, seq='CGTA')
print(dna.reverse_complement()) # DNASequence(len=4, seq='GCAT')

rna = dna.transcribe()
print(rna)                 # RNASequence('AUGC)

dna = DNASequence("ATGAAAATA")   # ATG AAA ATA -> M K I
protein = AminoAcidSequence.translate_aa(dna)
print(protein.seq)               # "MKI"
```
