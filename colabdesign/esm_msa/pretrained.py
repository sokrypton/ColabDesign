from colabdesign import esm_msa
import joblib


def get_model():
  alphabet = esm_msa.Alphabet.from_architecture('msa_transformer')
  config = esm_msa.config.CONFIG
  model = esm_msa.RunModel(alphabet, config)

  return model, alphabet
