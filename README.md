# CoMut
CoMut is a Python library for creating comutation plots to visualize genomic information.

<p align="center">
<img align="center" src="examples/images/melanoma_comut.svg" width="800"/>
</p>


## Installation

CoMut can be installed via `pip`:

`pip install comut`

## Quickstart

There is a [Quickstart notebook](https://github.com/vanallenlab/comut/blob/package/examples/quickstart.ipynb) that can create comut plots quickly with only input data. Simply download the notebook, specify the path to your data, and run the cells to produce a basic comut.

## Documentation

There is also a [Documentation notebook](https://github.com/vanallenlab/comut/blob/package/examples/documentation.ipynb) that documentation for CoMut. It describes the fundamentals of creating comuts, and provides the code used to generate the comut above.

## Dependencies

CoMut runs on python 3.7 or later. CoMut requires the following packages as dependencies (they will be installed along with CoMut if using `pip`)

```
numpy>=1.18.1
pandas>=0.25.3
palettable>=3.3.0
matplotlib>=3.3.1
```

## Versions

0.0.2 - Introduce compatability for Python 3.6
0.0.1 - Initial release