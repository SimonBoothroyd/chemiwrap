# Cheminformatics Wrappers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A proof of concept refactoring of the OpenFF toolkit registries and wrappers. This repo is not currently 
expected to be used by anyone, and is mainly created as a strawman.

See the [usage example to get started](examples/usage.py).

## Installation

The package and its dependencies can be installed using the `conda` package manager:

```shell
conda env create --name chemiwrap --file devtools/conda-ens/test-env.yaml
conda activate chemiwrap

python setup.py develop
```

If you have access to the OpenEye toolkits, it is recommended to install these also:

```shell
conda install -c openeye openeye-toolkits
```

#### License

This package is provided under the MIT license. It is derived from:

 * OpenFF Toolkit: Copyright (c) 2016-2019 Open Force Field Initiative

and is for the most part a refactored copy of the code found in the 
`openff.toolkits.utils.toolkits` module of the `openff-toolkit` package.

#### Copyright

Copyright (c) 2021, Simon Boothroyd