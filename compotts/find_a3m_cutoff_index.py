import csv
import numpy as np
import pandas as pd
from kneebow.rotor import Rotor
import pathlib

def find_a3m_cutoff_index(hhr_file, prettier_hhr_file_output=None):

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
        try:
            float(splitrow[2])
        except ValueError:
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

    # Sorting E-values
    df = pd.read_csv(prettier_hhr_file_output)
    df_evalues = df.sort_values(by=['E-value'])

    # Using kneebow to find knee
    y = list(df_evalues['E-value'])
    x = list(range(len(y)))
    data = np.array([[xi,yi] for xi,yi in zip(x,y)])

    rotor = Rotor()
    rotor.fit_rotate(data)
    elbow_idx = rotor.get_elbow_index()
    elbow_no = list(df_evalues['No'])[elbow_idx]


    return elbow_no
