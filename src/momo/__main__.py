from toml import load
from momo.sabdab import make_sabdab_dataset
from momo.skempi import make_skempi_dataset
from momo.config import MOMOConfig


def main():
    config_data = load("conf/config.toml")
    cfg = MOMOConfig(**config_data)

    skempi = make_skempi_dataset(
        filepath=cfg.filepaths.input.skempi,
        email=cfg.secret.email,
        ncbi_api_key=cfg.secret.ncbi_api_key,
    )
    print(skempi)
    sabdab = make_sabdab_dataset(
        filepath=cfg.filepaths.input.sabdab,
        email=cfg.secret.email,
        ncbi_api_key=cfg.secret.ncbi_api_key,
    )
    print(sabdab)


if __name__ == "__main__":
    main()
