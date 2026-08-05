import pandas as pd
import sys

df = pd.read_excel(sys.argv[1])
df = df.dropna(subset=['Transcript', 'HGVSc'])
df = df[df['HGVSc'].str.strip() != '']
hgvs = df['Transcript'] + ':' + df['HGVSc']
hgvs.to_csv(sys.argv[2], index=False, header=False)