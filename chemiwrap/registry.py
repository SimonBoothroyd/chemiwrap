from collections import defaultdict
from typing import Dict, List, Optional, Type, TypeVar, Union

from chemiwrap.providers import (
    AromaticityProvider,
    BondOrderProvider,
    ChargeProvider,
    ConformerProvider,
    FeatureProvider,
    FileIOProvider,
    SMARTSProvider,
    SMILESProvider,
    StereochemistryProvider,
    StereoisomerProvider,
    TautomerProvider,
)
from chemiwrap.providers.plugins import registered_feature_providers

FeatureProviderTypeVar = TypeVar("FeatureProviderTypeVar", bound=Type[FeatureProvider])

ConformerProviders = List[Union[str, Type[ConformerProvider]]]

ChargeProviders = List[Union[str, Type[ChargeProvider]]]
BondOrderProviders = List[Union[str, Type[BondOrderProvider]]]

AromaticityProviders = List[Union[str, Type[AromaticityProvider]]]
StereochemistryProviders = List[Union[str, Type[StereochemistryProvider]]]

TautomerProviders = List[Union[str, Type[TautomerProvider]]]
StereoisomerProviders = List[Union[str, Type[StereoisomerProvider]]]

SMARTSProviders = List[Union[str, Type[SMARTSProvider]]]
SMILESProviders = List[Union[str, Type[SMILESProvider]]]

FileIOProviders = List[Union[str, Type[FileIOProvider]]]

_REGISTRY_STACK = []
_DEFAULT_REGISTRY = None


def set_default_registry(feature_registry: "FeatureRegistry"):
    """Replace the default feature registry that should be used when no others have
    been activated using a context manager.

    Notes
    -----
    * This method can only be called when no other feature registries have been activated
      using a context manager.
    """

    if not isinstance(feature_registry, FeatureRegistry):
        raise TypeError("``feature_registry`` must be a ``FeatureRegistry`` instance")

    if len(_REGISTRY_STACK) > 0:

        raise RuntimeError(
            "The default feature can only be changed when it is also the current "
            "registry."
        )

    global _DEFAULT_REGISTRY
    _DEFAULT_REGISTRY = feature_registry


def current_registry() -> Optional["FeatureRegistry"]:
    """Returns the current feature registry that should be used."""
    return _DEFAULT_REGISTRY if len(_REGISTRY_STACK) == 0 else _REGISTRY_STACK[-1]


class FeatureRegistry:
    """A registry of 'features providers' that are used to offer the main functionality
    of the framework, such as conformer generation, partial charge generation or SMILES
    parsing.
    """

    @property
    def providers(self) -> Dict[str, List[Type[FeatureProvider]]]:
        """The current providers provided by this registry."""
        return {
            provider_type: [*provider_set]
            for provider_type, provider_set in self.__providers.items()
        }

    def __init__(
        self,
        conformer_providers: Optional[ConformerProviders] = None,
        charge_providers: Optional[ChargeProviders] = None,
        bond_order_providers: Optional[BondOrderProviders] = None,
        aromaticity_providers: Optional[AromaticityProviders] = None,
        stereochemistry_providers: Optional[StereochemistryProviders] = None,
        tautomer_providers: Optional[TautomerProviders] = None,
        stereoisomer_providers: Optional[StereoisomerProviders] = None,
        smarts_providers: Optional[SMARTSProviders] = None,
        smiles_providers: Optional[SMILESProviders] = None,
        file_io_providers: Optional[FileIOProviders] = None,
    ):
        """
        Parameters
        ----------
        conformer_providers
            The classes to use to provide conformer generation functionality in order
            of precedence. If ``None``, the providers will be inherited from the
            ``current_registry``.
            currently registered options.
        charge_providers
            The classes to use to provide partial charge generation functionality in
            order of precedence. If ``None``, the providers will be inherited from the
            ``current_registry``.
        bond_order_providers
            The classes to use to provide fractional bond order generation functionality
            in order of precedence. If ``None``, the providers will be inherited from the
            ``current_registry``.
        aromaticity_providers
            The classes to use to provide aromaticity models in order of precedence. If
            ``None``, the providers will be inherited from the ``current_registry``.
        stereochemistry_providers
            The classes to use to provide stereochemistry models in order of precedence.
             If ``None``, the providers will be inherited from the ``current_registry``.
        tautomer_providers
            The classes to use to provide tautomer enumeration functionality in order of
            precedence.  If ``None``, the providers will be inherited from the
            ``current_registry``.
        stereoisomer_providers
            The classes to use to provide stereoisomer enumeration functionality in
            order of precedence. If ``None``, the providers will be inherited from the
            ``current_registry``.
        smarts_providers
            The classes to use to provide SMARTS matching functionality in order of
            precedence. If ``None``, the providers will be inherited from the
            ``current_registry``.
        smiles_providers
            The classes to use to provide SMILES parsing / serialization functionality in
            order of precedence. If ``None``, the providers will be inherited from the
            ``current_registry``.
        file_io_providers
            The classes to use to provide file IO functionality in order of precedence.
            If ``None``, the providers will be inherited from the ``current_registry``.
        """

        kwargs = {
            # Ignore this hack for now...
            argument_name: provider_set
            for argument_name, provider_set in locals().items()
            if argument_name != "self"
        }

        self.__providers: Dict[str, List[Type[FeatureProvider]]] = {}

        providers_by_backend = {
            (provider_type, provider.backend()): provider
            for provider_type, providers in registered_feature_providers().items()
            for provider in providers
        }

        for argument_name, provider_set in kwargs.items():

            if provider_set is None:
                continue

            if len(provider_set) == 0:
                continue

            provider_type = "".join(x.capitalize() for x in argument_name.split("_"))[
                :-1
            ]

            provider_set = [
                providers_by_backend[(provider_type, provider)]
                if isinstance(provider, str)
                else provider
                for provider in provider_set
            ]

            self.__providers[provider_type] = provider_set

        if current_registry() is None:
            return

        for provider_type, provider_set in current_registry().providers.items():

            if provider_type in self.__providers:
                continue

            self.__providers[provider_type] = [*provider_set]

    def conformer_provider(self) -> Optional[ConformerProvider]:
        """Returns the provider that should be used to generate any conformers, or
        ``None`` if no such provider is available."""
        return self.get_provider_with_features(ConformerProvider, None)

    def charge_provider(
        self, model: Optional[Union[str, List[str]]] = None
    ) -> Optional[ChargeProvider]:
        """Returns the provider that should be used to generate partial charges for
        molecules, or ``None`` if no such provider is available.

        Parameters
        ----------
        model
            The (optional) charge model(s) that the providers should support.
        """

        return self.get_provider_with_features(
            ChargeProvider, [model] if isinstance(model, str) else model
        )

    def bond_order_provider(
        self, model: Optional[Union[str, List[str]]] = None
    ) -> Optional[BondOrderProvider]:
        """Returns the provider that should be used to generate fractional bond orders
        for molecules, or ``None`` if no such provider is available.

        Parameters
        ----------
        model
            The (optional) bond order model(s) that the providers should support.
        """

        return self.get_provider_with_features(
            BondOrderProvider, [model] if isinstance(model, str) else model
        )

    def aromaticity_providers(
        self, model: Optional[Union[str, List[str]]] = None
    ) -> Optional[AromaticityProvider]:
        """Returns the provider that should be used to determine the aromaticity
        of molecules, or ``None`` if no such provider is available.

        Parameters
        ----------
        model
            The (optional) aromaticity model(s) that the providers should support.
        """

        return self.get_provider_with_features(
            AromaticityProvider, [model] if isinstance(model, str) else model
        )

    def stereochemistry_provider(self) -> Optional[StereochemistryProvider]:
        """Returns the provider that should be used to determine the stereochemistry
        of molecules, or ``None`` if no such provider is available.
        """
        return self.get_provider_with_features(StereochemistryProvider, None)

    def tautomer_provider(self) -> Optional[TautomerProvider]:
        """Returns the provider that should be used to enumerate the tautomers of a
        molecule, or ``None`` if no such provider is available."""
        return self.get_provider_with_features(TautomerProvider, None)

    def stereoisomer_provider(self) -> Optional[StereoisomerProvider]:
        """Returns the provider that should be used to enumerate the stereoisomers of a
        molecule, or ``None`` if no such provider is available."""
        return self.get_provider_with_features(StereoisomerProvider, None)

    def smarts_provider(self) -> Optional[SMARTSProvider]:
        """Returns the provider that should be used to match SMARTS patterns against
        molecules, or ``None`` if no such provider is available."""
        return self.get_provider_with_features(SMARTSProvider, None)

    def smiles_provider(self) -> Optional[SMILESProvider]:
        """Returns the provider that should be used to interconvert molecules between
        SMILES representations, or ``None`` if no such provider is available."""
        return self.get_provider_with_features(SMILESProvider, None)

    def file_io_provider(self) -> Optional[FileIOProvider]:
        """Returns the provider that should be used to enumerate the tautomers of a
        molecule, or ``None`` if no such provider is available."""
        return self.get_provider_with_features(FileIOProvider, None)

    def get_provider_with_features(
        self, provider_type: FeatureProviderTypeVar, features: Optional[List[str]]
    ) -> Optional[FeatureProviderTypeVar]:

        if provider_type.__name__ not in self.__providers:
            return None

        providers = [
            provider
            for provider in self.__providers.get(provider_type.__name__, [])
            if provider.is_available()
            and (True if features is None else provider.supports_features(features))
        ]

        return None if len(providers) == 0 else providers[0]

    def get_provider_with_backend(
        self, provider_type: FeatureProviderTypeVar, backend: str
    ) -> Optional[FeatureProviderTypeVar]:

        if provider_type.__name__ not in self.__providers:
            return None

        for provider in self.__providers[provider_type.__name__]:

            if provider.backend().lower() != backend.lower():
                continue

            return provider

        return None

    @classmethod
    def using_backends(cls, *backends: str) -> "FeatureRegistry":
        """Creates a feature registry that contains providers which are based on the
        specified backends.

        Providers will be prioritised based on where their backend appears in the
        specified ``backends``.
        """

        backends = [backends] if isinstance(backends, str) else backends

        providers = defaultdict(list)

        for provider_type, provider_set in registered_feature_providers().items():

            provider_by_backend: Dict[str, Type[FeatureProvider]] = {
                provider.backend(): provider for provider in provider_set
            }

            for backend in backends:

                if backend not in provider_by_backend:
                    continue

                if not provider_by_backend[backend].is_available():
                    continue

                providers[provider_type].append(provider_by_backend[backend])

        registry = FeatureRegistry.__new__(FeatureRegistry)
        registry.__providers = {**providers}

        return registry

    def __enter__(self):
        _REGISTRY_STACK.append(self)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):

        if len(_REGISTRY_STACK) > 0:
            del _REGISTRY_STACK[-1]

    def __str__(self):
        def format_providers(providers: List[Type[FeatureProvider]]) -> str:
            return ", ".join(provider.backend() for provider in providers)

        return " ".join(
            f"{provider_type}=[{format_providers(providers)}]"
            for provider_type, providers in self.__providers.items()
        )

    def __repr__(self):
        return f"<FeatureSet {self.__str__()}>"


# TODO: replace with plugin-system + ordering by ~/.config/chemiwrap/config.yaml
set_default_registry(FeatureRegistry.using_backends("openeye", "rdkit", "ambertools"))
