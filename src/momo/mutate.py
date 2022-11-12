from copy import deepcopy
from dataclasses import dataclass, field
from typing import Callable

from pandas import DataFrame, Series


@dataclass
class Wildtype:
    pdb_code: str
    chains: dict[str, str] = field(repr=False)


@dataclass
class Mutation:
    old_residue: str
    chain: str
    position: int
    new_residue: str


@dataclass
class Variant:
    wildtype: Wildtype
    mutations: list[Mutation]
    id_: int
    chains: dict[str:str] = field(init=False)

    def __post_init__(self):
        self.chains = apply_mutations(self.wildtype.chains, self.mutations)


def apply_mutations(chains: dict[str:str], mutations: list[Mutation]) -> dict:
    mutated_chains = deepcopy(chains)
    for mutation in mutations:
        chain_key = next(key for key in chains.keys() if mutation.chain in key)
        mutated_chains[chain_key] = apply_substitution(
            chains[chain_key], mutation.position, mutation.new_residue
        )
    return mutated_chains


def apply_substitution(sequence: str, position: int, new_character: str) -> str:
    sequence_list = list(sequence)
    sequence_list[position - 1] = new_character
    return "".join(sequence_list)


def make_mutation(mutation_str: str) -> Mutation:
    return Mutation(
        old_residue=mutation_str[0],
        chain=mutation_str[1],
        position=int(mutation_str[2:-1]),
        new_residue=mutation_str[-1],
    )


def make_wildtype(pdb_code: str, chain_names: list[str], chains: list[str]) -> Wildtype:
    return Wildtype(pdb_code=pdb_code, chains=dict(zip(chain_names, chains)))


def make_variant(wildtype: Wildtype, mutations: list[Mutation], id_: int) -> Variant:
    mutation_list = [make_mutation(mutation) for mutation in mutations]
    return Variant(wildtype=wildtype, mutations=mutation_list, id_=id_)


def parse_mutation(df: DataFrame, parser: Callable) -> Series:
    return df["mutations"].apply(
        lambda mutations: [parser(mutation) for mutation in mutations]
    )


def variants_to_dataframe(variants: list[Variant]) -> DataFrame:
    variant_list = []
    for variant in variants:
        for chain_name, chain in variant.chains.items():
            mutations = {
                "old_residue": [],
                "chain": [],
                "position": [],
                "new_residue": [],
            }
            for mutation in variant.mutations:
                if mutation.chain in chain_name:
                    mutations["old_residue"].append(mutation.old_residue)
                    mutations["chain"].append(mutation.chain)
                    mutations["position"].append(mutation.position)
                    mutations["new_residue"].append(mutation.new_residue)
            variant_dict = {
                "pdb_code": variant.wildtype.pdb_code,
                "variant_id": variant.id_,
                "chain_name": chain_name,
                "mutated": any(
                    mutation.chain in chain_name for mutation in variant.mutations
                ),
                "variant_chain": chain,
                "wildtype_chain": variant.wildtype.chains[chain_name],
            } | mutations
            variant_list.append(variant_dict)
    return DataFrame(variant_list)
