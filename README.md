# Bioinformatic-supertool

В данном репозитории собраны инструменты для работы с нуклеотидными последовательностями.

**Точкой входа** в программу является файл `main.py`, в котором реализованы вызовы основных функций.

Скрипт `bio_files_processor.py` в котором реализованы функции `convert_multiline_fasta_to_oneline`, `parse_blast_output`.

- В папке `data` файл с fasta последовательностями (`example_multiline_fasta.fasta`), файл с последовательностями формата fastq: название, последовательность, "+" строка, строка качества (`example_fastq.fastq`), файл с результатом бласта (`example_blast_result.txt`), файл (`example_gbk.gbk`).
- В папке `utils` находятся модули с дополнительными функциями (`module_*.py`).
- В папке `filtered` хранятся результаты функций `filter_fastq` (**filt**), `convert_multiline_fasta_to_oneline` (**one_line.fasta**), `parse_blast_output` (**proteins.txt**).

## Функция `dna_rna_tools`

Позволяет выполнять базовые операции с последовательности ДНК/РНК:
- проверять является ли строка нуклеотидной последовательность;
- возвращать транскрибированную, обратную, комплементарную, обратную комплементарную последовательность.

Функция **dna_rna_tools** принимает на вход произвольное количество последовательностей нуклеотидов (`str`) и тип процесса, который надо произвести над последовательностями (`str`), 
`process` — тип операции (`"is_nucleic_acid"`, `"transcribe"`, `"reverse"`, `"complement"`, `"reverse_complement"`).

**Пример использования**

```python
run_dna_rna_tools('ATGC', 'is_nucleic_acid') # True
run_dna_rna_tools('AGU', 'CCuu' 'transcribe') # ['AGT', 'CCtt']
run_dna_rna_tools('ATg', 'reverse') # 'gTA'
run_dna_rna_tools('ctA', 'complement') # 'gaT'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
```

## Функция `filter_fastq`

Работает с последовательностями формата **fastq**. Позволяет фильтровать их по GC-составу, длине прочтения, качеству рида по шкале Phred33.

Функция **filter_fastq** принимает на вход 5 аргументов:
1. название файла с последовательностями формата fastq. 
2. название файла куда будет записан результат.
3. интервал (кортеж из двух значений) или значение (`float`) GC состава (в процентах) для фильтрации. В случае одной границы фильтруются (сохраняются) все риды ниже этой границы.
4. интервал (кортеж из двух значений) или значение (`float`) длины рида для фильтрации. В случае одной границы фильтруются все риды ниже этой границы.
5. пороговое значение среднего качества рида для фильтрации (`int`). Риды с качеством ниже порогового значения отбрасываются.

Сохраняет строки fastq-сиквенсов, отфильтрованных по заданным параметрам.

```python
filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0)
```

**Пример использования**
```python
filter_fastq('input_fastq', 'output_fastq', gc_bounds = (0,100),        length_bounds = (0,2**32), quality_threshold = 0)
```

## Функция `convert_multiline_fasta_to_oneline`

Читает поданный на вход fasta-файл, в котором последовательность (ДНК/РНК/белка) может быть разбита на несколько строк, после чего сохраняет в новый fasta-файл в котором каждая последовательность умещается в одну строку. 

Принимает на вход 2 аргумента (input_fasta и output_fasta). output_fasta необязательный аргумент, если не указан, возвращает список.

```python
convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str ='') -> list[str]
```

##  Функция `parse_blast_output`

Программа читает заданный txt файл, для каждого запроса QUERY (абзац Sequences producing significant alignments:), выбирает первую строку из столбца Description. Набор полученных белков сохраняет в алфавитном порядке в новый файл.

Принимает на вход 2 аргумента (input_file, output_file). 

```python
parse_blast_output(input_file: str, output_file: str) -> None
```
