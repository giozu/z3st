# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

from setuptools import find_packages, setup

setup(
    name="z3st",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[],
    python_requires=">=3.10",
    description="Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis",
    author="Giovanni Zullo",
    author_email="giovanni.zullo@polimi.it",
    license="Apache 2.0",
    url="https://github.com/giozu/z3st",
)
