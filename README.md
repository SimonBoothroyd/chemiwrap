# Cheminformatics Wrappers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A proof of concept refactoring of the OpenFF toolkit registries and wrappers. This repo is not currently 
expected to be used by anyone, and is mainly created as a strawman.

See the [usage example to get started](examples/usage.py).

### Installation

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

### Design Principles

The re-write of the toolkit wrappers was carried out with the following goals in mind:

* the 'features' provided by a toolkit (e.g. conformer generation, partial charge generation, etc.) should be exposed 
  in a uniform, strongly typed, and IDE friendly (including code completion and any static linting the **IDE** exposes)
  manner.
  
* new 'features' (e.g. rendering molecules) should be easy to implement **in the framework** in a way that is 
  independent of the way that others features are implemented.
  
* 'providers' of features (e.g. a class that exposed the ability to compute conformers using OpenEye) should be 
  implementable and discoverable by a plug-in system (see [plugins.py](chemiwrap/providers/plugins.py)).
  
* there should be minimum code duplication between features. An example of this is the validation of inputs to features
  (e.g. the maximum number of conformers to generate should be > 0) which is often independent of the backend providing 
  a feature.
  
* the provided 'features' that a 'feature provider' supports should be clear and programmatically queryable (e.g. what 
  charge models are supported by an OpenEye charge generator).

* anything produced by a 'feature provider' **must** have complete provenance attached. 

#### Feature Providers

At the core of this re-write is the notion of a *feature* and a *feature provider*.

* a *feature* is a specific high level functionality that the framework should provide. Examples of *features* include
  'generating conformers', 'computing partial charges', 'applying an aromaticity model'.
  
* a *feature provider* is a **class** in the framework that provides a features functionality

The framework provides an **abstract** base [``FeatureProvider`` class](chemiwrap/providers/_providers.py) that all 
feature providers **must** subclass. It provides a base interface that will be exposed by features including:

* ``backend`` - the backend that is providing a particular feature, e.g. "openeye"
* ``is_available`` - whether the feature provider is currently available e.g. can `rdkit` be imported.
* ``supports_features`` - a method that returns whether a particular feature (e.g. the AM1 charge model) is supported by 
  the provider.
  
The framework then also provides *10* additional **abstract base classes** that subclass ``FeatureProvider``:

* ``AromaticityProvider``
* ``BondOrderProvider``
* ``ChargeProvider``
* ``ConformerProvider``
* ``FileIOProvider``
* ``SMARTSProvider``
* ``SMILESProvider``
* ``StereochemistryProvider``
* ``StereoisomerProvider``
* ``TautomerProvider``

These **abstract base classes** represent the specific **features** that the framework supports, and define the 
**interface** for how the functionality they expose should be accessed. The classes themselves do not actually provide 
any functionality but rather act as templates for classes that will, *as well as implementing common duplicate 
functionality* such as input validation, attachment of default provenance, and more.  

Classes that actually implement a specific *feature* using a particular backend (e.g generating conformers using RDKit)
**must** subclass the relevant *feature provider* class, e.g. a class that provides functionality for generating 
conformers using OpenEye may look something like:

```python
class OpenEyeConformerProvider(ConformerProvider):
    """A wrapper around OpenEye's OMEGA conformer generation tool."""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    ...
```

Because `ConformerProvider` is flagged as *abstract* using the `abc.ABC` base class, your IDE / linter tell you
exactly what methods must be implemented in order for the class to be valid and usable. 

In the case of a `ConformerProvider` the main function that must be implemented is `_generate_conformers`

```python
@classmethod
def _generate_conformers(
    cls,
    molecule: Molecule,
    max_conformers: int,
    rms_tolerance: unit.Quantity = 1.0 * unit.angstrom,
) -> AnnotatedList[unit.Quantity]:

    ...
```

`_generate_conformers` is the *internal* implementation of the public `ConformerProvider.generate_conformers` function. 

The `ConformerProvider.generate_conformers` ([defined in `_providers.py`](chemiwrap/providers/_providers.py)) function 
handles validating the inputs (i.e. `max_conformers > 0`, `rms_tolerance > 0.0 Å` and has units compatible with 
angstroms), and then calls the `_generate_conformers` function implemented by the subclass (using built-in the 
`classmethod` machinery) to handle the actual conformer generation. See for example the [actual `OpenEyeConformerProvider`
implementation](chemiwrap/providers/openeye.py), which can me used according to:

```python
from openff.toolkit.topology import Molecule
from chemiwrap.providers.openeye import OpenEyeConformerProvider

conformers = OpenEyeConformerProvider.generate_conformers(
    molecule=Molecule.from_smiles("C"), max_conformers=10
)
# More on this later...
print(conformers.provenance)
```

Splitting the main functions exposed by a *feature provider* into a public and internal function allows the base class 
to handle all the logic that would be duplicated across multiple implementations, thus making the actual 
implementation a lot simpler and easier to maintain.  

Using class inheritance gives for free:

* docstring inheritance without the need for additional decorators / copy mechanisms.
* a uniform, well typed, consistent API across all providers of a feature that can be 
  easily linted by the commonly used IDEs as well as tools like `mypy`.

#### Feature Registry

A key design on the OpenFF toolkit is the `ToolkitRegistry` class. The toolkit registry is responsible for storing 
multiple toolkit wrappers and providing a single API point for calling the functions provided by each contained wrapper,
such that the registry attempts to call the requested function on each contained wrapper in a pre-specified order 
until either a wrapper returns a non-error return value or there are no more wrappers to call on.

Arguably this 'waterfall' of function calls is a limitation of how toolkit wrappers are designed. Because a toolkit 
registry does not know what functionality each function of a toolkit wrapper supports (e.g. what charge models can be 
applied), it must try each wrapper it knows about until it finds one that does support the requested functionality.

The knock-on effect of this is that if say a toolkit wrapper yields a genuine error (i.e. not an error saying this 
functionality is not supported) when applying its functionality (e.g. the OpenEye wrapper raises an exception because 
the molecule passed to assign partial charges is garbage), the registry will assume that the wrapper simply does not 
support the functionality and move onto the next one. Depending on the implementation of the next wrapper, it may allow 
the 'garbage molecule' that caused the problem to be used through without exception, and ultimately return a 'garbage 
output'.

The [``FeatureRegistry`` class](chemiwrap/providers/_providers.py) aims to provide a similar functionality while avoiding
the problem outlined here. It is a container of ``FeatureProvider`` classes. It can currently **only** contain the 10 
feature providers outlined above.

It can be initialized either using specific implementations of the different feature providers at a 'per feature' level:

```python
from chemiwrap.registry import FeatureRegistry

feature_registry = FeatureRegistry(
    conformer_providers=[OpenEyeConformerProvider, RDKitConformerProvider],
    charge_providers=[AmberToolsChargeProvider],
    ...
)
```

or by specifying the `backend` that should provide individual features:

```python
feature_registry = FeatureRegistry(
    conformer_providers=["openeye", "rdkit"],
    charge_providers=["ambertools"]
)
```

or by specifying the overall `backend`'s that should provide the features that they support: 

```python
feature_registry = FeatureRegistry.using_backends("openeye", "rdkit")
```

in all cases the order in which the feature providers (/ backends) are specified is the order of priority in which each 
*feature provider* should be used.

The feature provider *with the highest priority* can be accessed directly from the registry:

```python
feature_registry.conformer_provider().generate_conformers(...)
```

The 'feature provider accessors' are methods rather than properties. This is because, for certain feature providers, the 
user can request the feature provider with the highest priority that supports specific functionality, e.g.,

```python
charge_provider = feature_registry.charge_provider(model=["am1bcc", "am1bccelf10"])
```

This is a critical design point. Because the user can access the provider that providers the **exact** functionality 
they need, a registry **does not** support 'try all of the available functions until you find one that doesn't 
raise an exception'. Arguably in *almost all* cases this is closer to the desired behaviour, and should lead to less 
'user suprise'. Namely, an exception is raised immediately when a function fails out at a 'bad' input, rather than 
propagating input down the call the chain and likely silently producing 'garbage output'.

Because all the returned providers are *strongly typed* and not *dynamically generated* everything works
seamlessly with IDEs (especially their built-in linters and **code completion**). In principle the feature registry 
class (especially the init signature) could be dynamically generated in the future when IDE support for dynamic 
classes is more robust and allows for better code completion / linting.

#### Global Feature Registries

Much like the OpenFF toolkit, the proposed re-write includes support for creating a 'global' feature registry that 
should be used by all function calls that require a particular feature.

This can be accessed by calling `current_registry()`, e.g.

```python
from chemiwrap.registry import current_registry

current_registry().conformer_provider().generate_conformers(...)
```

The initial default registry is currently hardcoded and depends on which providers are currently available, but in 
future will be configurable using a `~/.config/chemiwrap/config.yaml` file (or a config pointed to by an env variable).
It can be programmatically changed by calling `set_default_registry`.

The re-write goes further than the OpenFF toolkit by also allowing the user to create a local scope in which a particular
feature registry should be used rather than the default global one. This is made possible using context managers:

```python
with FeatureRegistry(
    conformer_providers=["openeye"],
    charge_providers=["openeye"],
) as feature_registry:

    assert current_registry() == feature_registry  # True

    # We can access features using the local context manager variable.
    conformers = (
        feature_registry()
        .conformer_provider()
        .generate_conformers(
            molecule, max_conformers=10, rms_tolerance=1.0 * unit.angstrom
        )
    )

    # OR access the the 'global' current registry which will be equivilant to
    # ``feature_registry`` in this scope.
    conformers = (
        current_registry()
        .conformer_provider()
        .generate_conformers(
            molecule, max_conformers=10, rms_tolerance=1.0 * unit.angstrom
        )
    )
```

If only a subset of providers are passed to a ``FeatureRegistry``, e.g.,

```python
FeatureRegistry(conformer_providers=[OpenEyeConformerProvider])
```

the remaining providers (e.g. ``charge_providers``) will be inherited from the ``current_registry()``.
If you don't want this to be the case, you will either need to specify all feature providers explicitly or 
pass ``[]`` to the providers that should not be made available.

Feature registries can be infinitely nested, see [for example `usage.py`](examples/usage.py)

In general, code that requires a particular *feature* should **always** access the associated provider through the 
``current_registry()`` API.

#### Provenance

A big change introduced in this re-write is that **any output** returned by a *feature provider* will
have provenance information attached. Currently, the return values are one of ``AnnotatedList``, ``AnnotatedDict``, or
``AnnotatedQuantity`` (see [here](chemiwrap/utilities/provenance.py)).

These 'annotated' classes behave exactly like lists, dicts, or quantities with the addition that they have 
a 'provenance' attribute.

An example of this is generating partial charges:

```python
molecule = Molecule.from_smiles("CCCC")

with FeatureRegistry(conformer_providers=["openeye"]):

    conformers = (
        current_registry()
        .conformer_provider()
        .generate_conformers(
            molecule, max_conformers=10, rms_tolerance=1.0 * unit.angstrom
        )
    )
    pprint(conformers.provenance)
```

which yields an output of 

```text
{'backend': 'openeye',
 'max-conformers': 10,
 'openeye-version': '2020.2.4',
 'rms-tolerance': '1.0 A'}
```

see [usage.py](examples/usage.py) for more examples of the provenance returned.

### License

This package is provided under the MIT license. It is derived from:

 * OpenFF Toolkit: Copyright (c) 2016-2019 Open Force Field Initiative

and is for the most part a refactored copy of the code found in the 
`openff.toolkits.utils.toolkits` module of the `openff-toolkit` package.

### Copyright

Copyright (c) 2021, Simon Boothroyd