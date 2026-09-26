from setuptools import setup, find_packages
setup(
    name='brain-pasta',
    version='0.0.0',
    python_requires=">=3.9",
    packages=find_packages(),
    install_requires=[
        'numpy>=1.26,<2.0',
        'scipy>=1.13',
        'scikit-learn>=1.6',
        'scikit-learn-extra>=0.3',
        'setuptools'
    ],
)