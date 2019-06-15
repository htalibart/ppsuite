import pandas as pd

def get_info(column_name, info_res_file):
    df = pd.read_csv(info_res_file)
    return df[column_name][0]
