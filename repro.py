# %%
import tiledb
import pandas as pd
import numpy as np
from os.path import exists
import shutil

# %%
ARRAY_NAME = "feature_counts"

if exists(ARRAY_NAME):
    shutil.rmtree(ARRAY_NAME, ignore_errors=True)

# create the dimensions sequence x gene
dims = [
    tiledb.Dim(name="sequence", tile=None, dtype="ascii"),
    tiledb.Dim(name="gene", tile=None, dtype="ascii"),
]

# create the domain
dom = tiledb.Domain(*dims)

# create the count attribute (set the default filler to 0)
count_attr = tiledb.Attr(name="count", dtype=np.int32, fill=0)

# create the array schema
schema = tiledb.ArraySchema(domain=dom, sparse=True, attrs=[count_attr])

tiledb.Array.create(ARRAY_NAME, schema)

# %%
#
# this cell just reads the data from counts (which might have one row, or might have many rows)
#
df = pd.read_csv("counts.tsv", sep="\t", header=0)
genes = np.array(df.columns[1:])
samples = np.array(df["sequencerunid"])

for sample in samples:
    counts = df.loc[df["sequencerunid"] == sample].values.flatten().tolist()[1:]
    sample_coords = np.repeat(sample, len(counts))

    print(len(counts))
    print(len(genes))
    print(len(sample))

    with tiledb.open(ARRAY_NAME, mode="w") as A:
        A[sample_coords, genes] = counts

# %%
#
# this cell compresses the above with inline data
#
g = ["A1BG", "A1BG-AS1", "A1CF", "A2M", "A2M-AS1"]
c = [14, 17, 994, 2596, 47]
s = np.repeat("01595556-0bb1-4416-8da9-cdca78fd4ab4", len(g))

print(len(c))
print(len(g))

with tiledb.open(ARRAY_NAME, mode="w") as A:
    A[s, g] = c
