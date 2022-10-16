import jax

# adopted from https://github.com/deepmind/alphafold/blob/main/alphafold/model/prng.py
class SafeKey:
  """Safety wrapper for PRNG keys."""

  def __init__(self, key):
    self._key = key
    self._used = False

  def _assert_not_used(self):
    if self._used:
      raise RuntimeError('Random key has been used previously.')

  def get(self):
    self._assert_not_used()
    self._used = True
    return self._key

  def split(self, num_keys=2):
    self._assert_not_used()
    self._used = True
    new_keys = jax.random.split(self._key, num_keys)
    return jax.tree_map(SafeKey, tuple(new_keys))

  def duplicate(self, num_keys=2):
    self._assert_not_used()
    self._used = True
    return tuple(SafeKey(self._key) for _ in range(num_keys))
