from setuptools import setup, find_packages

setup(
    name="codon_usage_analysis",
    version="0.1.0",
    description="A package for analyzing codon usage and calculating CAI values",
    author="Ruixuan",
    packages=find_packages(),
    install_requires=[
        "biopython",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)