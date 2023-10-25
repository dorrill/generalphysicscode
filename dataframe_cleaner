# Code for Reading in and Working with Pandas Data Frames, Cleaning up the Data, and then analyzing it
import numpy as np
import pandas as pd
pd.__version__
#Url of data frame: example, housing dataset
url = 'https://github.com/mattharrison/datasets/raw/master/data/ames-housing-dataset.zip'
url = 'data/ames-housing-dataset.zip'
df = pd.read_csv(url, engine='pyarrow', dtype_backend='pyarrow')

# Data cleaning functions
def shrink_ints(df):
    mapping = {}
    for col in df.dtypes[df.dtypes=='int64[pyarrow]'].index:
        max_ = df[col].max()
        min_ = df[col].min()
        if min_ < 0:
            continue
        if max_ < 255:
            mapping[col] = 'uint8[pyarrow]'
        elif max_ < 65_535:
            mapping[col] = 'uint16[pyarrow]'
        elif max_ <  4294967295:
            mapping[col] = 'uint32[pyarrow]'
    return df.astype(mapping)


def clean_data(df):
    return (df
     .assign(**df.select_dtypes('string').replace('', 'Missing').astype('category'),
             **{'Garage Yr Blt': df['Garage Yr Blt'].clip(upper=df['Year Built'].max())})
     .pipe(shrink_ints)
    )    

clean_DF = clean_data(df).dtypes
