from dataclasses import dataclass, field
from pandas import DataFrame, Series
from typing import Callable


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
    chains: dict[str, str] = field(init=False)

    def __post_init__(self):
        self.chains = apply_mutations(self.wildtype.chains, self.mutations)


def apply_mutations(chains: dict, mutations: list[Mutation]) -> dict:
    mutated_chains = {}
    for mutation in mutations:
        chain_key = next(key for key in chains.keys() if mutation.chain in key)
        mutated_chains[chain_key] = apply_substitution(
            chains[chain_key], mutation.position, mutation.new_residue
        )
    return mutated_chains


def apply_substitution(sequence: str, position: int, new_character: str) -> str:
    sequence_list = list(sequence)
    sequence_list[position] = new_character
    return "".join(sequence_list)


def make_mutation(mutation_str: str) -> Mutation:
    return Mutation(
        old_residue=mutation_str[0],
        chain=mutation_str[1],
        position=int(mutation_str[2:-1]),
        new_residue=mutation_str[-1],
    )


def make_wildtype(
    pdb_code: str, chain_keys: list[str], sequences: list[str]
) -> Wildtype:
    return Wildtype(pdb_code=pdb_code, chains=dict(zip(chain_keys, sequences)))


def make_variant(wildtype: Wildtype, mutations: list[Mutation], id_: int) -> Variant:
    mutation_list = [make_mutation(mutation) for mutation in mutations]
    return Variant(wildtype=wildtype, mutations=mutation_list, id_=id_)


# def mutate(sequence: str, positions: list[int], new_characters: list[str], row) -> str:
#     sequence_list = list(sequence)
#     for position, new_character in zip(positions, new_characters):
#         sequence_list[int(position) - 1] = new_character
#     return "".join(sequence_list)


# def apply_mutation(df: DataFrame) -> DataFrame:
#     return df.assign(
#         variant=df.apply(
#             lambda x: mutate(x["sequence"], x["position"], x["new_residue"])
#             if x["mutation_chain"].isin(x["chain"])
#             else x["sequence"],
#             axis=1,
#         )
#     )


# def map_mutations(df: DataFrame) -> DataFrame:
#     return df.pipe(split_mutations).pipe(assign_mutations).pipe(apply_mutation)
