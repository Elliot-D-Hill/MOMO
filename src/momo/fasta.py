from momo.utils import collapse_on_column


class FastaParser:
    def __init__(self):
        self.column_names = [
            "pdb_code",
            "chain",
            "description",
            "sequence",
        ]

    def process_header(self, df):
        header = df["header"].str.split("|")
        df["pdb_code"] = header.str[1]
        description = header.str[2]
        df["chain"] = description.str[0]
        df["description"] = description.str.split(", ").str[-1]
        return df.drop("header", axis=1)

    def organize_dataframe(self, df):
        return df[self.column_names].sort_values(
            ["pdb_code", "chain"], ignore_index=True
        )

    def run_pipeline(self, df):
        return (
            df.pipe(self.process_header)
            .pipe(collapse_on_column, "chain")
            .pipe(self.organize_dataframe)
        )
