import importlib
import sys
from collections import defaultdict
from typing import Dict, List, Type

from chemiwrap.providers import FeatureProvider, __feature_types__

if sys.version_info < (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points

_RESERVED_BACKENDS = ["ambertools", "rdkit", "openeye", "built-in"]


def registered_feature_providers() -> Dict[str, List[Type[FeatureProvider]]]:

    # Import the built-in providers so they register with __subclass__
    for module_name in ["ambertools", "openeye", "rdkit", "builtin"]:
        importlib.import_module(f"chemiwrap.providers.{module_name}")

    for entry_point_modules in entry_points(group="chemiwrap.providers"):

        for entry_point_module in entry_point_modules:
            entry_point_module.load()

    feature_providers: Dict[str, List[Type[FeatureProvider]]] = {
        feature_type.__name__: [*feature_type.__subclasses__()]
        for feature_type in __feature_types__
    }

    # Perform validation to make sure the built-in backends aren't overwritten.
    for provider_type, providers in feature_providers.items():

        backend_providers = defaultdict(list)

        for provider in providers:
            backend_providers[provider.backend()].append(provider)

        for backend in _RESERVED_BACKENDS:

            external_providers = [
                provider
                for provider in backend_providers[backend]
                if not provider.__module__.startswith("chemiwrap.providers")
            ]

            if len(external_providers) == 0:
                continue

            raise RuntimeError(
                f"The external {external_providers[0].__module__} module registers a "
                f"`{provider_type}` feature provider with a {backend} backend. "
                f"This backend is currently reserved to prevent accidentally "
                f"overwriting the built-in feature providers. Consider modifying the "
                f"backend of the provider to 'custom-{backend}'."
            )

    return feature_providers
