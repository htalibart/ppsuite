import csv
import numpy as np
import pandas as pd
from kneebow.rotor import Rotor
import pathlib
from Bio.Blast import NCBIXML


def get_elbow_no(df, column_name='E-value'):
    # Sorting E-values
    df_evalues = df.sort_values(by=[column_name])

    # Using kneebow to find knee
    y = list(df_evalues[column_name])
    x = list(range(len(y)))
    data = np.array([[xi,yi] for xi,yi in zip(x,y)])

    rotor = Rotor()
    rotor.fit_rotate(data)
    elbow_idx = rotor.get_elbow_index()
    elbow_no = list(df_evalues['No'])[elbow_idx]
    return elbow_no


def list_can_be_cast_to_float(l):
    can = True
    for x in l:
        try:
            float(x)
        except ValueError:
            can = False
    return can


def find_hhblits_cutoff_index(hhr_file, prettier_hhr_file_output=None):

    if prettier_hhr_file_output is None:
        hhr_file_name = str(hhr_file)
        prettier_hhr_file_output = pathlib.Path('.'.join(hhr_file_name.split('.')[:-1])+"_hhr.csv")

    # Reading .hhr file
    with open(hhr_file, 'r') as f:
        lines = f.readlines()

    startline=1000000000000000
    stopline=1000000000000000
    for nline in range(len(lines)):
        if lines[nline].startswith("  1"):
            startline = nline
        if lines[nline]=="\n" and nline>startline and nline<stopline:
            stopline = nline

    rows = []
    interesting_lines = lines[startline:stopline]
    for line in interesting_lines:
        splitrow = line.split(' ')
        splitrow = list(filter(None, splitrow))
        while not list_can_be_cast_to_float(splitrow[2:-3]):
            splitrow[1]+=splitrow[2]
            del splitrow[2]
        row = splitrow[:6]
        rows.append(row)

    # Turning it into a parsable .csv file
    header = "No,Hit,Prob,E-value,P-value,Score\n"
    with open(prettier_hhr_file_output, 'w') as f:
        f.write(header)
        csvwriter = csv.writer(f)
        for row in rows:
            csvwriter.writerow(row)

    df = pd.read_csv(prettier_hhr_file_output)
    return get_elbow_no(df, column_name='E-value')



def find_blast_cutoff_index(blast_results_file, blast_csv=None):

    if blast_csv is None:
        blast_results_file_name = str(blast_results_file)
        blast_csv = pathlib.Path('.'.join(blast_results_file_name.split('.')[:-1])+'.csv')
    
    rows = []
    no = 1
    with open(blast_results_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for recordcount, blast_record in enumerate(blast_records):
                for alignment, desc in zip(blast_record.alignments,blast_record.descriptions):
                    for hsp in alignment.hsps:
                        hit = desc.title.split(' ')[0]
                        row = [no, hit, hsp.expect]
                        rows.append(row)
                        no+=1

    header = "No,Hit,E-value\n"
    with open(blast_csv, 'w') as f:
        f.write(header)
        csvwriter = csv.writer(f)
        for row in rows:
            csvwriter.writerow(row)

    df = pd.read_csv(blast_csv)
    return get_elbow_no(df, column_name='E-value')
