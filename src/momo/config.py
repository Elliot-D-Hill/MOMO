from typing import Optional
from pydantic import BaseModel
from pathlib import Path


class Paths(BaseModel):
    intermediate: Path
    input: Path
    output: Path
    models: Path


class InputFiles(BaseModel):
    skempi: Path
    sabdab: Path


class IntermediateFiles(BaseModel):
    skempi_fasta: Path
    sabdab_fasta: Path


class OutputFiles(BaseModel):
    skempi: Path
    sabdab: Path
    chain_classifier: Path


class Files(BaseModel):
    input: InputFiles
    intermediate: IntermediateFiles
    output: OutputFiles


class TrainingParameters(BaseModel):
    learning_rate: float
    batch_size: int
    n_epochs: int
    train_proportion: float
    clip_grad: float
    lr_scheduler_gamma: float


class ModelParameters(BaseModel):
    embedding_size: int


class DataParameters(BaseModel):
    vocabulary_size: int
    num_classes: int


class Parameters(BaseModel):
    training: TrainingParameters
    model: ModelParameters
    data: DataParameters


class Secret(BaseModel):
    email: str
    ncbi_api_key: str


class MOMOConfig(BaseModel):
    files: Files
    paths: Paths
    filepaths: Optional[Files]
    parameters: Parameters
    secret: Secret

    def __init__(self, **data) -> None:
        super().__init__(**data)
        self.make_filepaths()

    def make_filepaths(self) -> None:
        input_files = InputFiles(
            skempi=self.paths.input / self.files.input.skempi,
            sabdab=self.paths.input / self.files.input.sabdab,
        )
        intermediate_files = IntermediateFiles(
            skempi_fasta=self.paths.intermediate / self.files.intermediate.skempi_fasta,
            sabdab_fasta=self.paths.intermediate / self.files.intermediate.sabdab_fasta,
        )
        output_files = OutputFiles(
            skempi=self.paths.output / self.files.output.skempi,
            sabdab=self.paths.output / self.files.output.sabdab,
            chain_classifier=self.paths.models / self.files.output.chain_classifier,
        )
        self.filepaths = Files(
            input=input_files, intermediate=intermediate_files, output=output_files
        )
