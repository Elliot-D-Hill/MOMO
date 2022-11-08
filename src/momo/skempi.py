from pandas import DataFrame
from momo.utils import clean_dataframe_header


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
            .pipe(self.organize_dataframe)
        )
