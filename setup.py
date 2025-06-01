from setuptools import setup, find_packages

setup(
    name='ngs_pipeline',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # dependency
    ],
    entry_points={
        'console_scripts': [
            'rna-seq = scripts.rna_seq.rna_seq:main',
	    'chip-seq = scripts.chip_seq.chip_seq:main'
        ],
    },
)
