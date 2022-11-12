from dataclasses import asdict
from toml import load
from pandas import DataFrame, read_csv

from proteintools.mutate import map_mutations
from proteintools.embedding_model import SequenceClassifier
from proteintools.fastatools import (
    make_dataframe_from_fasta,
    get_fasta_from_ncbi_query,
    FastaParser,
)
from momo.mutate import variants_to_dataframe
from momo.sabdab import SabdabParser
from momo.skempi import SkempiParser, make_variants, make_wildtypes
from momo.config import MOMOConfig


def main():
    config_data = load("conf/config.toml")
    cfg = MOMOConfig(**config_data)

    # SKEMPI
    skempi_df = read_csv(cfg.filepaths.input.skempi, sep="\t")
    skempi_parser = SkempiParser()
    skempi_parsed = skempi_parser.run_pipeline(skempi_df)
    print(skempi_parsed)
    skempi_pdb_codes = skempi_parsed["pdb_code"].unique()
    skempi_fasta_text = get_fasta_from_ncbi_query(
        skempi_pdb_codes, cfg.secret.email, cfg.secret.ncbi_api_key
    )
    skempi_fasta_df = make_dataframe_from_fasta(skempi_fasta_text)
    fasta_parser = FastaParser()
    skempi_fasta_parsed = fasta_parser.run_pipeline(skempi_fasta_df)
    print(skempi_fasta_parsed)

    wildtypes = make_wildtypes(skempi_fasta_parsed)
    wildtypes_dict = {wt.pdb_code: wt for wt in wildtypes}
    variants = make_variants(skempi_parsed, wildtypes_dict)
    variant_df = variants_to_dataframe(variants)
    skempi = variant_df.merge(skempi_parsed, on=["pdb_code", "variant_id"], how="inner")
    print(skempi)

    # SABDAB
    sabdab_df = read_csv(cfg.filepaths.input.sabdab, sep="\t")
    sabdab_pdb_codes = sabdab_df["pdb"].unique()
    sabdab_fasta_text = get_fasta_from_ncbi_query(
        sabdab_pdb_codes, cfg.secret.email, cfg.secret.ncbi_api_key
    )
    sabdab_fasta_df = make_dataframe_from_fasta(sabdab_fasta_text)
    print(sabdab_fasta_df)
    sabdab_fasta_parsed = fasta_parser.run_pipeline(sabdab_fasta_df)
    sabdab_parser = SabdabParser()
    sabdab_parsed = sabdab_parser.run_pipeline(sabdab_df, sabdab_fasta_parsed)
    print(sabdab_parsed)


if __name__ == "__main__":
    main()
