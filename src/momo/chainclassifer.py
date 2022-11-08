from toml import load

from pandas import read_csv
from torch import device, cuda
from torch.nn import CrossEntropyLoss
from torch.optim import AdamW, lr_scheduler

from proteintools.embedding_model import (
    EmbeddingModel,
    SequenceClassifier,
    build_vocabulary,
)
from momo.config import MOMOConfig


def main():
    config_data = load("conf/config.toml")
    cfg = MOMOConfig(**config_data)

    # FIXME deep to reference nested parameters
    params = cfg.parameters
    df = read_csv(cfg.filepaths.output.sabdab)
    X = df["sequence"]
    y = df["chain_type"]
    print(y.shape, X.shape)
    # switch device to GPU, if it is available
    processing_device = device("cuda" if cuda.is_available() else "cpu")
    # initialize model objects
    vocabulary = build_vocabulary(X)
    vocab_size = len(vocabulary)
    num_class = y.nunique()
    model = SequenceClassifier(vocab_size, params.embedding_size, num_class).to(
        processing_device
    )
    criterion = CrossEntropyLoss()
    optimizer = AdamW(model.parameters(), lr=params.learning_rate)
    scheduler = lr_scheduler.StepLR(optimizer, 1.0, gamma=params.lr_scheduler_gamma)
    chain_classifier = EmbeddingModel(
        X,
        y,
        vocabulary,
        params,
        model,
        criterion,
        optimizer,
        scheduler,
        processing_device,
    )
    # fit model
    test_accuracy = chain_classifier.run_pipeline()
    chain_classifier.save_model(cfg.filepaths.chain_classifier)
    n_train = int(len(y) * chain_classifier.train_proportion)
    n_test_val = int(len(y) * ((1 - chain_classifier.train_proportion) / 2))
    # print info
    print("test accuracy {:8.3f}".format(test_accuracy))
    print(f"Trained on {n_train} sequences")
    print(f"Validated on {n_test_val} sequences")
    print(f"Tested on {n_test_val} sequences")


if __name__ == "__main__":
    main()
