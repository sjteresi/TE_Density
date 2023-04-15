from setuptools import setup, find_packages


setup(
    name="te-density",
    version="2.1.1",
    description="Calculates Transposable Element density",
    url="https://github.com/sjteresi/TE_Density",
    packages=find_packages(),
    author="Scott Teresi, Michael Teresi",
    license = "GPL-3.0",
    install_requires=[
        'coloredlogs>=15.0',
        'h5py>=3.7',
        'matplotlib>=3.6',
        'numpy>=1.23',
        'numexpr>=2.8.3',
        'pandas>=1.5',
        'scipy>=1.9',
        'tqdm>=4.64',
    ],
    scripts=[
        './process_genome.py'
    ],
)
