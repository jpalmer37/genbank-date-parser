# Genbank Date Parser
Simple Python date parser. Intended for phylogenetic tree dating, especially with LSD2 in IQTree. 


```
python gb_date_parser.py --help
usage: gb_date_parser.py [-h] [-d HEADER_DELIM] [-p ACCNO_POS] fasta

positional arguments:
  fasta                 Input FASTA file for which dates should be retrieved.

optional arguments:
  -h, --help            show this help message and exit
  -d HEADER_DELIM, --header-delim HEADER_DELIM
                        Delimiter character used in the FASTA header.
  -p ACCNO_POS, --accno-pos ACCNO_POS
                        Zero-indexed position of the accession number in the FASTA header
```
