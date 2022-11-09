from numpy import select


class SabdabParser:
    def __init__(self):
        self.column_names = [
            "pdb_code",
            "chain",
            "chain_type",
            "description",
            "sequence",
        ]

    def format_dataframe(self, df):
        old_columns = ["pdb", "Hchain", "Lchain", "antigen_chain"]
        new_columns = ["pdb_code", "heavy_chain", "light_chain", "antigen_chain"]
        columns = dict(zip(old_columns, new_columns))
        return (
            df.filter(items=old_columns)
            .assign(pdb=df["pdb"].str.upper())
            .assign(
                antigen_chain=df["antigen_chain"].str.replace(" | ", "", regex=False)
            )
            .rename(columns=columns)
        )

    def filter_dataframe(self, df):
        return df[df["heavy_chain"] != df["light_chain"]]

    def group_sums(self, df):
        return df.groupby("pdb_code").sum()

    def merge_dataframes(self, df, fasta_df):
        return df.merge(fasta_df, on="pdb_code", how="inner")

    def assign_chain_type(self, df):
        choices = ["heavy_chain", "light_chain", "antigen_chain"]
        conditions = [is_chain_subset(df, chain) for chain in choices]
        df["chain_type"] = select(conditions, choices, default="antigen_chain")
        return df

    def run_pipeline(self, sabdab_df, fasta_df):
        return (
            sabdab_df.pipe(self.format_dataframe)
            .pipe(self.filter_dataframe)
            .pipe(self.group_sums)
            .pipe(self.merge_dataframes, fasta_df)
            .pipe(self.assign_chain_type)
        )


# TODO refactor this and make more abstract; "chain" is too specific
def is_chain_subset(df, column):
    return [any([i in str(a) for i in str(b)]) for a, b in zip(df[column], df["chain"])]
