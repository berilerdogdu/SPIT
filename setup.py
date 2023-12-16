from setuptools import setup, find_packages

setup(
    name='spit',
    version='0.1.8',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'matplotlib>=3.4.3',
        'matplotlib_venn>=0.11.6',
        'numpy>=1.20.3',
        'pandas>=1.3.4',
        'rpy2>=3.5.14',
        'scikit_learn>=1.1.1',
        'scipy>=1.11.1',
        'seaborn>=0.11.2',
        'tqdm>=4.62.3',
    ],
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'spit = spit.run_spit:main',  # Entry point for script1
        ],
    },
    package_data={'mypackage': ['r_scripts/*.R']},
)
