import pandas as pd
import numpy as np
import sys
from scipy.stats import chi2


fname = sys.argv[1]


df = pd.read_csv(fname, sep="\t")
print(df.head())
df["P"] = chi2.sf(np.square(df.Z), 1)

df.to_csv("{fname}.fixed".format(fname=fname), index=False, sep="\t")
