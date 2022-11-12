from pandas import DataFrame, Series, read_csv
from momo.mutate import make_variant, make_wildtype, variants_to_dataframe
from momo.utils import clean_dataframe_header
from proteintools.fastatools import (
    make_dataframe_from_fasta,
    get_fasta_from_ncbi_query,
    FastaParser,
)


class SkempiParser:
    def __init__(self) -> None:
        self.old_names = [
            "pdb",
            "mutations_cleaned",
            "affinity_wt_parsed",
            "affinity_mut_parsed",
        ]
        self.new_names = ["pdb_code", "mutation", "wildtype_kd", "variant_kd"]
        self.header_replacements = {" ": "_", "-": "_", "(": "", ")": "", "#": ""}

    def format_dataframe(self, df: DataFrame) -> DataFrame:
        columns = dict(zip(self.old_names, self.new_names))
        return clean_dataframe_header(df, self.header_replacements).rename(
            columns=columns
        )

    def filter_rows(self, df: DataFrame) -> DataFrame:
        antibodies = ["AB/AG", "AB/AG,Pr/PI"]
        is_antibody = df["hold_out_type"].isin(antibodies)
        return df[is_antibody]

    def filter_columns(self, df: DataFrame) -> DataFrame:
        return df[self.new_names]

    def trim_pdb_code(self, df: DataFrame) -> DataFrame:
        return df.assign(pdb_code=df["pdb_code"].str[0:4]).reset_index(drop=True)

    def assign_variant_id(self, df: DataFrame) -> DataFrame:
        df["variant_id"] = df.groupby("pdb_code").cumcount().add(1)
        return df

    def split_mutations(self, df: DataFrame) -> Series:
        df["mutation"] = df["mutation"].str.split(",")
        return df

    # FIXME .dropna() here?
    def organize_dataframe(self, df):
        self.new_names.insert(1, "variant_id")
        return df[self.new_names].sort_values(["pdb_code", "variant_id"])

    def run_pipeline(self, df: DataFrame) -> DataFrame:
        return (
            df.pipe(self.format_dataframe)
            .pipe(self.filter_rows)
            .pipe(self.filter_columns)
            .pipe(self.trim_pdb_code)
            .pipe(self.assign_variant_id)
            .pipe(self.split_mutations)
            .pipe(self.organize_dataframe)
        )


def make_wildtypes(df):
    return df.groupby("pdb_code").apply(
        lambda group: make_wildtype(group.name, group["chain"], group["sequence"])
    )


def make_variants(df, wildtypes_dict):
    return df.apply(
        lambda x: make_variant(
            wildtype=wildtypes_dict[x["pdb_code"]],
            mutations=x["mutation"],
            id_=x["variant_id"],
        ),
        axis=1,
    )


def make_skempi_dataset(filepath, email, ncbi_api_key):
    skempi_df = read_csv(filepath, sep="\t")
    skempi_parser = SkempiParser()
    skempi_parsed = skempi_parser.run_pipeline(skempi_df)
    skempi_pdb_codes = skempi_parsed["pdb_code"].unique()
    skempi_fasta_text = get_fasta_from_ncbi_query(skempi_pdb_codes, email, ncbi_api_key)
    skempi_fasta_df = make_dataframe_from_fasta(skempi_fasta_text)
    fasta_parser = FastaParser()
    skempi_fasta_parsed = fasta_parser.run_pipeline(skempi_fasta_df)
    wildtypes = make_wildtypes(skempi_fasta_parsed)
    wildtypes_dict = {wt.pdb_code: wt for wt in wildtypes}
    variants = make_variants(skempi_parsed, wildtypes_dict)
    variant_df = variants_to_dataframe(variants)
    skempi = variant_df.merge(skempi_parsed, on=["pdb_code", "variant_id"], how="inner")
    return skempi
