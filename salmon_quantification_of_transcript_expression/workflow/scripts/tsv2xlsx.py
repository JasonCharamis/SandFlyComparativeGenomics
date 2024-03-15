import csv
import re
import argparse
from xlsxwriter.workbook import Workbook

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

tsv_file = args.filename
xlsx_file = re.sub('.tsv',".xlsx",tsv_file)
  
workbook = Workbook(xlsx_file)
worksheet = workbook.add_worksheet()

read_tsv = csv.reader(open(tsv_file, 'r'), delimiter='\t')
  
for row, data in enumerate(read_tsv):
    worksheet.write_row(row, 0, data)
  
workbook.close()
