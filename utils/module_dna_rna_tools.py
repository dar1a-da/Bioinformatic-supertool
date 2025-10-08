def check_nucleic_acid(seq: str) -> bool:
    return set(seq).issubset('aAtTcCgG') or set(seq).issubset('aAuUcCgG')


def transcribe_nucleic_acid(seq: str) -> str:
    trans_rules = {'t': 'u', 'T': 'U', 'u': 't', 'U': 'T'}
    return ''.join(trans_rules.get(el, el) for el in seq)


def reverse_nucleic_acid(seq: str) -> str:
    return seq[::-1]


def complement_nucleic_acid(seq: str) -> str:
    compl_rules = {'a': 't', 'A': 'T', 't': 'a', 'T': 'A',
                   'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}
    return ''.join(compl_rules.get(el, el) for el in seq)


def reverse_complement_nucleic_acid(seq: str) -> str:
    return reverse_nucleic_acid(complement_nucleic_acid(seq))