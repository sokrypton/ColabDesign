from setuptools import setup, find_packages
setup(
    name='ColabDesign',
    version='1.0.1',
    packages=find_packages(include=['colabdesign',
                                    'colabdesign.af','colabdesign.tr',
                                    'colabdesign.alphafold']),
    install_requires=['py3Dmol','absl-py','biopython',
                      'chex','dm-haiku','dm-tree',
                      'immutabledict','jax','ml-collections',
                      'numpy','pandas','scipy','tensorflow'],
)
