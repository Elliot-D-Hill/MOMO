from pandas import DataFrame


def clean_dataframe_header(df: DataFrame, replacements: dict) -> DataFrame:
    df.columns = df.columns.str.strip().str.lower()
    df.rename({})
    for key, value in replacements.items():
        df.columns = df.columns.str.replace(key, value, regex=False)
    return df


def collapse_on_column(df: DataFrame, column_name: str) -> DataFrame:
    """Collapse rows that differ only by specified column"""
    columns = [name for name in df.columns if name != column_name]
    return (
        df.groupby(columns)[column_name]
        .apply(lambda lst: "".join(sorted(lst)))
        .reset_index()
    )
