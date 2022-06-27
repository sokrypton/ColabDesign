from setuptools import setup, find_packages
setup(
    name='ColabDesign',
    version='1.0.1',
    packages=find_packages(include=['colabdesign','colabdesign.af','colabdesign.tr']),
    install_requires=['py3Dmol'],
)
