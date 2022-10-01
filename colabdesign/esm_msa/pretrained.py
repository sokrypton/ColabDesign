import esm_jax
import joblib


def get_model():
    alphabet = esm_jax.Alphabet.from_architecture('msa_transformer')
    config = esm_jax.config.CONFIG
    model = esm_jax.RunModel(alphabet, config)

    return model, alphabet
