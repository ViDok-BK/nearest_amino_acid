import numpy as np
import pandas as pd

path = "C:\\Users\\Admin\\OneDrive\\Desktop\\Github\\result.csv"

df = pd.read_csv(path)

frame = df.groupby(["residue_name", "residue_number"])
with open("result.txt", 'w') as f:
    for item in frame:
        f.write(item[0])