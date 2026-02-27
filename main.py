from abc import ABC, abstractmethod
from typing import ClassVar, Mapping

class InvalidAlphabetError(ValueError):
    pass

class BiologicalSequence(ABC):
    """
    Work with the len function.
    The ability to get elements by index and make sequence slices.
    Printing in a beautiful way.
    The ability to check the alphabet of the sequence for correctness.
    """
    ALPHABET: ClassVar[frozenset[str]]

    def __init__(self, seq: str):
        self.seq = seq
        self.check_alphabet()

    def __len__(self) -> int:
        return len(self.seq)
    
    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.__class__(self.seq[item])
        return self.seq[item]
    
    def __iter__(self):
        return iter(self.seq)
    
    def __str__(self) -> str:
        cls = self.__class__.__name__
        MAX_LEN = 60
        if len(self.seq) <= MAX_LEN:
            preview = self.seq
        else:
            preview = self.seq[:MAX_LEN - 3] + "..."
        return f"{cls}(len={len(self)}, seq='{preview}')"
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(seq='{self.seq}')"
        

    @abstractmethod
    def check_alphabet(self) -> None:
        raise NotImplementedError
    

class NucleicAcidSequence(BiologicalSequence):
    
    COMPLEMENT_MAP: ClassVar[Mapping[str, str]]
    ALPHABET: ClassVar[frozenset[str]]

    def check_alphabet(self):
        bad = set(self.seq) - self.ALPHABET
        if bad:
            raise InvalidAlphabetError(
                f"{self.__class__.__name__}: недопустимые символы {bad}")

    def complement(self) -> 'NucleicAcidSequence':
        cmap = self.__class__.COMPLEMENT_MAP
        return self.__class__(''.join(cmap.get(ch, ch) for ch in self.seq))
    
    def reverse(self) -> 'NucleicAcidSequence':
        return self.__class__(self.seq[::-1])

    def reverse_complement(self) -> 'NucleicAcidSequence':
        return self.complement().reverse()
    

class DNASequence(NucleicAcidSequence):
    
    ALPHABET: ClassVar[frozenset[str]] = frozenset(set('aAtTgGcC'))
    COMPLEMENT_MAP: ClassVar[Mapping[str, str]] = {'a': 't', 'A': 'T',
                                                   't': 'a', 'T': 'A',
                                                   'g': 'c', 'G': 'C',
                                                   'c': 'g', 'C': 'G'}

    def transcribe(self) -> 'RNASequence':
        trans_rules = {'t': 'u', 'T': 'U'}
        rna = ''.join(trans_rules.get(ch, ch) for ch in self.seq)
        return RNASequence(rna)


class RNASequence(NucleicAcidSequence):
   
    ALPHABET: ClassVar[frozenset[str]] = frozenset(set('aAuUgGcC'))
    COMPLEMENT_MAP: ClassVar[Mapping[str, str]] = {'a': 'u', 'A': 'U',
                                                   'u': 'a', 'U': 'A',
                                                   'g': 'c', 'G': 'C',
                                                   'c': 'g', 'C': 'G'}


class AminoAcidSequence(BiologicalSequence):
    
    ALPHABET: ClassVar[frozenset[str]] = frozenset(
        set('ACDEFGHIKLMNPQRSTVWY*Xacdefghiklmnpqrstvwy*x'))
    
    CODON_TABLE: ClassVar[dict[str, str]] = {
        # Phenylalanine
        "UUU": "F", "UUC": "F",
        # Leucine
        "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        # Isoleucine
        "AUU": "I", "AUC": "I", "AUA": "I",
        # Methionine (start)
        "AUG": "M",
        # Valine
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        # Serine
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
        # Proline
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        # Threonine
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        # Alanine
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        # Tyrosine
        "UAU": "Y", "UAC": "Y",
        # Histidine
        "CAU": "H", "CAC": "H",
        # Glutamine
        "CAA": "Q", "CAG": "Q",
        # Asparagine
        "AAU": "N", "AAC": "N",
        # Lysine
        "AAA": "K", "AAG": "K",
        # Aspartic acid
        "GAU": "D", "GAC": "D",
        # Glutamic acid
        "GAA": "E", "GAG": "E",
        # Cysteine
        "UGU": "C", "UGC": "C",
        # Tryptophan
        "UGG": "W",
        # Arginine
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        # Glycine
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        # Stop codons
        "UAA": "*", "UAG": "*", "UGA": "*",}
    
    def check_alphabet(self):
        bad = set(self.seq) - self.ALPHABET
        if bad:
            raise InvalidAlphabetError(f"{self.__class__.__name__}: недопустимые символы {bad}")
        
    @classmethod
    def translate_aa(cls, nuc: 'NucleicAcidSequence',
                     frame: int=0,
                     stop_at_stop: bool=True,
                     unknown: str='X') -> 'AminoAcidSequence':
        """Converts DNA/RNA into amino acids by triplets.

        Params:
        nuc : NucleicAcidSequence
            DNASequence or RNASequence.
        frame : int
            Frame shift (0, 1 or 2).
        stop_at_stop : bool
            If True, we stop at the first stop codon (*).
            If False, we insert '*' and continue.
        unknown : str
            How to replace codons with ambiguous characters/unknown codons.
        Returns:
        AminoAcidSequence
        """
        if frame not in (0, 1, 2):
            raise ValueError('frame must be 0, 1 or 2')
        
        seq = nuc.seq.upper()
        seq = seq.replace('T', 'U')

        protein = []

        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]

            aa = cls.CODON_TABLE.get(codon, unknown)

            if aa == '*' and stop_at_stop:
                break

            protein.append(aa)

        return cls(''.join(protein))

