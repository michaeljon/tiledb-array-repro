# %%
import tiledb
import pandas as pd
import numpy as np

from tqdm import tqdm

import shutil
from os.path import exists

# Print TileDB-Py version
print("TileDB-Py version: " + tiledb.__version__)

# Print TileDB core version
print("TileDB core version: " + str(tiledb.libtiledb.version()))

ARRAY_NAME = "feature_counts"

# %%
if exists(ARRAY_NAME):
    shutil.rmtree(ARRAY_NAME, ignore_errors=True)

# create the dimensions sequence x gene
dims = [
    tiledb.Dim(name="sequencerunid", tile=None, domain=("", ""), dtype="ascii"),
    tiledb.Dim(name="gene", tile=None, domain=("", ""), dtype="ascii"),
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
df = pd.read_csv("counts-with-sequence-rows.tsv.gz", sep="\t", header=0, index_col=0, compression="gzip")

genes = df.columns
samples = df.index

with tiledb.open(ARRAY_NAME, mode="w") as A:
    for ndx, sample in enumerate(tqdm(samples)):
        counts = df.iloc[ndx, :]
        sample_coords = np.repeat(sample, len(counts))

        A[sample_coords, genes] = counts

# %%
with tiledb.open(ARRAY_NAME, mode="r") as A:
    df = A.multi_index[
        [
            "01595556-0bb1-4416-8da9-cdca78fd4ab4",
            "0265cfd4-d874-49cf-b53e-61e7998547c4",
            "02d7e850-1fc6-4c3a-8779-f073ab910180",
            "03350543-6ed3-43be-a4fd-a08ac1472c8f",
            "03f8bcb2-a096-4a99-b389-e09d3803ca5b",
            "0486d86e-077b-41e3-ac41-f8df6d8c209f",
            "04b3a558-418f-4f3b-abeb-859531486dc0",
            "055d56b5-8422-4671-bf86-555a80c25ee0",
            "061bbd2b-3732-4220-9253-9099a0b7ad48",
        ],
        [
            "A2M",
            "A2M-AS1",
            "A2ML1",
        ],
    ]

    print(df["count"].mean())

    # print(A.df[A.df["gene"].str.startswith("MIR")].groupby("gene").mean("count"))
# %%
