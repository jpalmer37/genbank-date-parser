# %%
import os, sys
import pandas as pd
from Bio import SeqIO, Entrez
from glob import glob
from datetime import datetime as dt
import argparse
import re

#%%
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', help="Input FASTA file for which dates should be retrieved.")
	# parser.add_argument('-i', '--inplace', action='store_true', help="Indicates whether to append the dates to the end of the FASTA header with a | delimiter")
	parser.add_argument('-d', '--header-delim', default="|", type=str, help="Delimiter character used in the FASTA header.")
	parser.add_argument('-p', '--accno-pos', default=0, type=int, help="Zero-indexed position of the accession number in the FASTA header")

	return parser.parse_args()
	

def parse_dates(accno_list):
		
	query = Entrez.efetch(db='nucleotide', id=accno_list, rettype='gb', retmode='text')

	collection_dates = {}

	expr_accno = re.compile('LOCUS\s+([^\s]+)\s+')
	expr_date = re.compile('=\"(.+)\"')

	# if collection dates are not present, use the oldest submitted journal date as closest estimate
	expr_date_backup = re.compile("Submitted \((.+)\)")

	accno = ''
	backup_date = ''
	for line in query:
		line = line.strip("\n ")

		if line.startswith("LOCUS"):
			# if the main collection date could not be found and the backup date is available,
			# then use the backup date instead 
			if accno != '' and backup_date != '' and collection_dates[accno] == '':
				collection_dates[accno] = backup_date

			# start of the next gb file 
			accno = expr_accno.search(line).group(1)
			collection_dates[accno] = ''
			backup_date = ''

		# collect and parse backup dates from journal entries
		if line.startswith('JOURNAL   Submitted'):
			date_search = expr_date_backup.search(line).group(1)
			date_search = dt.strptime(date_search, '%d-%b-%Y')

			if isinstance(backup_date, str):
				backup_date = date_search
			else:
				if date_search < backup_date:
					backup_date = date_search

		# retrieve collection dates (preferred date choice)
		if line.startswith("/collection_date"):
			date = expr_date.search(line).group(1)
			collection_dates[accno] = date
		
	# types of genbank date formats
	# written this way with extendability in mind
	# (expression, parse_format, output_format)
	expr_dates = [
		(re.compile('\d{2}-[A-z]+-\d{4}'), '%d-%b-%Y', '%Y-%m-%d'),
		(re.compile('^[A-z]+-\d{4}'), '%b-%Y', '%Y-%m'),
		(re.compile('^\d{4}$'), '%Y', '%Y'),
	]

	for accno, date in collection_dates.items():
		if isinstance(date, dt):
			collection_dates[accno] = date.strftime('%Y-%m-%d')

		elif isinstance(date, str):
			for expr, in_fmt, out_fmt in expr_dates:
				if expr.match(date):
					collection_dates[accno] = dt.strptime(date, in_fmt).strftime(out_fmt)
					break  
		else:
			print(accno)
			print(date)
			print("NOT VALID FORMAT")

	return collection_dates

#%%
def main():
	args = init_parser()

	# parse sequences in FASTA format 
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	# retrieve the accession numbers from the FASTA file  
	accno_list = [x.id.split(args.header_delim)[args.accno_pos] for x in seqs]

	# retrieve the collection dates 
	dates_dict = parse_dates(accno_list)

	output_basename = os.path.basename(args.fasta).split(".")[0]

	# print the dates TSV file 
	date_df = pd.DataFrame(list(dates_dict.items()), columns=['accno','date'])
	date_df.to_csv(output_basename + "_dates.tsv", sep='\t', index=False)

if __name__ == '__main__':
	main()