from pandas import DataFrame


def clean_dataframe_header(df: DataFrame, replacements: dict) -> DataFrame:
    df.columns = df.columns.str.strip().str.lower()
    df.rename({})
    for key, value in replacements.items():
        df.columns = df.columns.str.replace(key, value, regex=False)
    return df


def save_data(df, path, filename, postfix):
    path.mkdir(parents=True, exist_ok=True)
    prefix = filename.split(".")[0]
    df.to_csv(f"{path}/{prefix}_{postfix}.csv", index=False)
