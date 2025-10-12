import os

def read_fastq(file_fastq: str) -> dict[str, tuple[str, str]]:
    """
    Reads and converts in dict (parse) fastq file

    :param file_fastq: Путь до файла с NGS разметкой
    :return seqs: Словарь секса
    """
    with open(file_fastq) as file:
        lines = file.readlines()

    seqs = {}  
    for i in range(0, len(lines), 4):
        key = lines[i].strip('\n')
        seq = lines[i+1].strip('\n')
        plus = lines[i+2].strip('\n')
        qual = lines[i+3].strip('\n')
        seqs[key] = (seq, plus, qual)

    return seqs

def write_fastq(seqs, file):
    os.makedirs('filtered', exist_ok=True)
    with open(os.path.join('filtered',file), mode='w') as file:
        for i in seqs.keys():
            file.write(i+'\n'+seqs[i][0]+'\n'+seqs[i][1]+'\n'+seqs[i][2]+'\n')
