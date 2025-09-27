from setuptools import setup, find_packages

setup(
    name='spit',
    version='2.0.2',
    description='Comprehensive Python package for detecting differential transcript usage (DTU) in RNA-seq data.',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib>=3.4.3',
        'matplotlib_venn>=0.11.6',
        'numpy>=1.20.3',
        'pandas>=1.3.4',
        'scikit_learn>=1.1.1',
        'scipy>=1.11.1',
        'seaborn>=0.11.2',
        'tqdm>=4.62.3',
    ],
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'spit = spit.run_spit:main',
        ],
    },
)
