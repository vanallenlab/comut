import setuptools

with open("README.md") as f:
    LONG_DESCRIPTION, LONG_DESC_TYPE = f.read(), "text/markdown"

with open("reqs/base-requirements.txt") as f:
    REQUIREMENTS = f.read().splitlines()

NAME = "comut"
AUTHOR_NAME, AUTHOR_EMAIL = "Jett Crowdis", "jcrowdis@broadinstitute.org"
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
LICENSE = "MIT"
DESCRIPTION = "A Python library for creating comutation plots to visualize genomic information"
URL = "https://github.com/vanallenlab/comut/tree/package"
PYTHON_REQ = ">=3.6"

setuptools.setup(
    name=NAME,
    version="0.0.2",
    author=AUTHOR_NAME,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESC_TYPE,
    url=URL,
    packages=setuptools.find_packages(),
    classifiers=CLASSIFIERS,
    python_requires=PYTHON_REQ,
    install_requires=REQUIREMENTS,
)
