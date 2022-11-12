from numpy import select
from pandas import DataFrame, read_csv
from proteintools.fastatools import (
    make_dataframe_from_fasta,
    get_fasta_from_ncbi_query,
    FastaParser,
)


class SabdabParser:
    def __init__(self):
        self.column_names = [
            "pdb_code",
            "chain",
            "chain_type",
            "description",
            "sequence",
        ]

    def assign_chain_type(self, df):
        choices = ["heavy_chain", "light_chain", "antigen_chain"]
        conditions = [is_chain_subset(df, chain) for chain in choices]
        df["chain_type"] = select(conditions, choices, default="antigen_chain")
        return df

    def run_pipeline(self, sabdab_df, fasta_df):
        old_columns = ["pdb", "Hchain", "Lchain", "antigen_chain"]
        new_columns = ["pdb_code", "heavy_chain", "light_chain", "antigen_chain"]
        columns = dict(zip(old_columns, new_columns))
        return (
            sabdab_df.filter(items=old_columns)
            .rename(columns=columns)
            .assign(
                pdb_code=lambda df: df["pdb_code"].str.upper(),
                antigen_chain=lambda df: df["antigen_chain"].str.replace(
                    " | ", "", regex=False
                ),
            )
            .pipe(lambda df: df[df["heavy_chain"] != df["light_chain"]])
            .groupby("pdb_code")
            .sum()
            .merge(fasta_df, on="pdb_code", how="inner")
            .pipe(self.assign_chain_type)
        )


# TODO refactor this and make more abstract; "chain" is too specific
def is_chain_subset(df, column):
    return [any([i in str(a) for i in str(b)]) for a, b in zip(df[column], df["chain"])]


def make_sabdab_dataset(filepath: str, email: str, ncbi_api_key: str) -> DataFrame:
    sabdab_df = read_csv(filepath, sep="\t")
    sabdab_pdb_codes = sabdab_df["pdb"].unique()
    sabdab_fasta_text = get_fasta_from_ncbi_query(sabdab_pdb_codes, email, ncbi_api_key)
    sabdab_fasta_df = make_dataframe_from_fasta(sabdab_fasta_text)
    fasta_parser = FastaParser()
    sabdab_fasta_parsed = fasta_parser.run_pipeline(sabdab_fasta_df)
    sabdab_parser = SabdabParser()
    sabdab = sabdab_parser.run_pipeline(sabdab_df, sabdab_fasta_parsed)
    return sabdab
