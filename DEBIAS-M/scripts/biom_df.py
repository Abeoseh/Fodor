import biom



def load_table(path: str, h: bool=False):
    """pandas df from biom table"""
    table = biom.load_table(path)

    df = table.to_dataframe()
    df = df.rename_axis("IDs").reset_index()
    if h:
        print(df.head(n=10))
    
    return df

load_table("./BIOM/80310/otu_table.biom")