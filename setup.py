from setuptools import setup, find_packages
setup(
    name='colabdesign',
    version='1.1.0',
    packages=find_packages(include=['colabdesign*']),
    install_requires=['py3Dmol','absl-py','biopython',
                      'chex','dm-haiku','dm-tree',
                      'immutabledict','jax','ml-collections',
                      'numpy','pandas','scipy','optax','joblib'],
    include_package_data=True
)
