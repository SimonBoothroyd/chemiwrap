import abc
import copy
from typing import Dict, List, Optional, Tuple, Union

from openff.toolkit.topology import Molecule
from simtk import unit
from typing_extensions import TypedDict

from chemiwrap.exceptions import (
    IncompatibleUnitsError,
    InvalidArgumentError,
    MissingFeatureProviderError,
    UnsupportedFeatureError,
)
from chemiwrap.utilities import deduplicate_molecules
from chemiwrap.utilities.provenance import (
    AnnotatedDict,
    AnnotatedList,
    AnnotatedQuantity,
)
from chemiwrap.utilities.typing import AtomStereochemistry, BondStereochemistry

DefaultConformerSettings = TypedDict(
    "DefaultConformerSettings", {"min_conformers": int, "default_conformers": int}
)


class FeatureProvider(abc.ABC):
    """A base for classes that provide a specific feature provided by the toolkit,
    e.g. conformer generation, charge assignment, etc."""

    @classmethod
    @abc.abstractmethod
    def backend(cls) -> str:
        """A unique name associated with the underlying backend for the features
        provided by this class, e.g. 'built-in', 'openeye', 'rdkit'."""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def is_available(cls) -> bool:
        """Returns whether the underlying backend is available, e.g. is the required
        backend module importable and licensed."""
        raise NotImplementedError()

    @classmethod
    def supports_features(cls, features: Union[str, List[str]]) -> bool:
        """Returns whether this provider supports a specific feature. This may be, for
        example, a specific charge or aromaticity model."""

        features = [features] if isinstance(features, str) else features

        if len(features) == 0:
            return True

        raise UnsupportedFeatureError(cls, ", ".join(features), None)


class ModelProvider(FeatureProvider, abc.ABC):
    """A base for classes that provide a specific model (e.g. a charge or aromaticity
    model) that can be applied to a molecule."""

    @classmethod
    @abc.abstractmethod
    def supported_models(cls) -> Tuple[str, ...]:
        """Returns a list of the models that can be applied by this class."""
        return tuple()

    @classmethod
    def supports_features(cls, features: Union[str, List[str]]) -> bool:
        """Returns whether this provider supports a specific feature. This may be, for
        example, a specific charge or aromaticity model."""

        features = [features] if isinstance(features, str) else features
        return all(feature in cls.supported_models() for feature in features)


class ConformerProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for generating conformers
    for molecules."""

    @classmethod
    def _generate_conformers(
        cls,
        molecule: Molecule,
        max_conformers: int,
        rms_tolerance: unit.Quantity = 1.0 * unit.angstrom,
    ) -> AnnotatedList[unit.Quantity]:
        """The internal implementation of the `generate_conformers` function.

        Notes
        -----
        * It is safe to assume that the inputs have been validated.
        """
        raise NotImplementedError()

    @classmethod
    def generate_conformers(
        cls,
        molecule: Molecule,
        max_conformers: int,
        rms_tolerance: unit.Quantity = 1.0 * unit.angstrom,
    ) -> AnnotatedList[unit.Quantity]:
        """Generates a set of conformers for the provided molecule.

        Parameters
        ----------
        molecule:
            The molecule to generate conformers for.
        max_conformers:
            The maximum number of conformers to generate.
        rms_tolerance:
            If the minimum RMSD between any generated conformers.

        Returns
        -------
            A list of the generated conformers, where each conformer is a unit wrapped
            numpy array with shape=(n_atoms, 3)
        """

        if not rms_tolerance.unit.is_compatible(unit.angstrom):

            raise IncompatibleUnitsError(
                "rms_tolerance", rms_tolerance.unit, unit.angstrom
            )

        if rms_tolerance < 0.0 * unit.angstrom:
            raise InvalidArgumentError("``rms_tolerance`` must be > 0 Å")

        if max_conformers <= 0:
            raise InvalidArgumentError("``max_conformers`` must be > 0")

        return cls._generate_conformers(
            copy.deepcopy(molecule), max_conformers, rms_tolerance
        )

    @classmethod
    def assign_conformers(
        cls,
        molecule: Molecule,
        max_conformers: int = 1,
        rms_tolerance: unit.Quantity = 1.0 * unit.angstrom,
    ):
        """Generates a set of conformers and assigns them to the provided molecule
        in-place.

        Parameters
        ----------
        molecule:
            The molecule to generate conformers for.
        max_conformers:
            The maximum number of conformers to generate.
        rms_tolerance:
            If the minimum RMSD between any generated conformers.
        """

        molecule._conformers = cls.generate_conformers(
            molecule, max_conformers, rms_tolerance
        )


class ChargeProvider(ModelProvider, abc.ABC):
    """The base for classes that provide functionality for computing partial charges for
    molecules.
    """

    @classmethod
    @abc.abstractmethod
    def supported_models(cls) -> Tuple[str, ...]:
        """Returns a list of the models that can be used by this class to compute
        partial charges."""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def default_conformer_settings(
        cls,
    ) -> Dict[str, Optional[DefaultConformerSettings]]:
        """The default settings to use when generating conformers when none have
        been provided for the computation.
        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def _generate_partial_charges(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Union[List[unit.Quantity], AnnotatedList[unit.Quantity]],
    ) -> AnnotatedQuantity:
        """The internal implementation of the `compute_partial_charges` function."""
        raise NotImplementedError()

    @classmethod
    def generate_partial_charges(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Optional[
            Union[List[unit.Quantity], AnnotatedList[unit.Quantity]]
        ] = None,
    ) -> AnnotatedQuantity:
        """Computes a set of partial charges for the specified molecule using the
        requested model.

        Parameters
        ----------
        molecule
            The molecule to compute partial charges for.
        model
            The charge model to use. See ``supported_models`` for
            details.
        conformers
            The (optional) set of conformers to average the partial charges over. If
            none are provided a set will automatically be generated. These should be
            unit wrapped numpy arrays with shape=(n_atoms, 3).

        Returns
        -------
            The set of computed partial charges.
        """
        from chemiwrap.registry import current_registry

        if model not in cls.supported_models():
            raise UnsupportedFeatureError(cls, model, cls.supported_models())

        if conformers is None:

            default_conformer_settings = cls.default_conformer_settings()[model]

            if default_conformer_settings is None:
                conformers = []
            else:

                conformer_provider = current_registry().conformer_provider()

                if conformer_provider is None:

                    raise MissingFeatureProviderError(
                        ConformerProvider.__name__,
                        None,
                        f"A conformer provider must be available when generating "
                        f"{model} charges if a list of `conformers` is not provided.",
                    )

                conformers = conformer_provider.generate_conformers(
                    molecule, default_conformer_settings["default_conformers"]
                )

        return cls._generate_partial_charges(copy.deepcopy(molecule), model, conformers)

    @classmethod
    def assign_partial_charges(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Optional[
            Union[List[unit.Quantity], AnnotatedList[unit.Quantity]]
        ] = None,
    ):
        """Computes a set of partial charges for the specified molecule using the
        requested model, and directly assigns them to the molecule in-place.

        Parameters
        ----------
        molecule
            The molecule to compute partial charges for.
        model
            The charge model to use. See ``supported_models`` for
            details.
        conformers
            The (optional) set of conformers to average the partial charges over. If
            none are provided a set will automatically be generated. These should be
            unit wrapped numpy arrays with shape=(n_atoms, 1).
        """

        molecule.partial_charges = cls.generate_partial_charges(
            molecule, model, conformers
        )


class BondOrderProvider(ModelProvider, abc.ABC):
    """The base for classes that provide functionality for computing fractional
    bond orders for molecules."""

    @classmethod
    @abc.abstractmethod
    def supported_models(cls) -> Tuple[str, ...]:
        """Returns a list of the models that can be used by this class to compute
        fractional bond orders."""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def default_conformer_settings(
        cls,
    ) -> Dict[str, Optional[DefaultConformerSettings]]:
        """The default settings to use when generating conformers when none have
        been provided for the computation.
        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def _generate_fractional_bond_orders(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Union[List[unit.Quantity], AnnotatedList[unit.Quantity]],
    ) -> AnnotatedList[float]:
        """The internal implementation of the ``generate_fractional_bond_orders``
        function.
        """
        raise NotImplementedError()

    @classmethod
    def generate_fractional_bond_orders(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Optional[
            Union[List[unit.Quantity], AnnotatedList[unit.Quantity]]
        ] = None,
    ) -> AnnotatedList[float]:
        """Computes a set of fractional bond orders for the specified molecule using the
        requested model.

        Parameters
        ----------
        molecule
            The molecule to compute fractional bond orders for.
        model
            The bond order model to use. See ``supported_models``
            for details.
        conformers
            The (optional) set of conformers to average the fractional bond orders over.
            If none are provided a set will automatically be generated. These should be
            unit wrapped numpy arrays with shape=(n_atoms, 3).

        Returns
        -------
            The set of computed fractional bond orders with shape=(n_atoms,).
        """

        from chemiwrap.registry import current_registry

        if model not in cls.supported_models():
            raise UnsupportedFeatureError(cls, model, cls.supported_models())

        if conformers is None:

            default_conformer_settings = cls.default_conformer_settings()[model]

            if default_conformer_settings is None:
                conformers = []
            else:

                conformer_provider = current_registry().conformer_provider()

                if conformer_provider is None:
                    raise MissingFeatureProviderError(
                        ConformerProvider.__name__,
                        None,
                        f"A conformer provider must be available when generating "
                        f"{model} fractional bond orders if a list of `conformers` is "
                        f"not provided.",
                    )

                conformers = conformer_provider.generate_conformers(
                    molecule, default_conformer_settings["default_conformers"]
                )

        return cls._generate_fractional_bond_orders(
            copy.deepcopy(molecule), model, conformers
        )

    @classmethod
    def assign_fractional_bond_orders(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Optional[
            Union[List[unit.Quantity], AnnotatedList[unit.Quantity]]
        ] = None,
    ):
        """Computes a set of fractional bond orders for the specified molecule using the
        requested model and directly assigns them to the molecule in-place.

        Parameters
        ----------
        molecule
            The molecule to compute fractional bond orders for.
        model
            The bond order model to use. See the ``bond_order_methods`` for details.
        conformers
            The (optional) set of conformers to average the fractional bond orders over.
            If none are provided a set will automatically be generated. These should be
            unit wrapped numpy arrays with shape=(n_atoms, 1).
        """

        bond_orders = cls.generate_fractional_bond_orders(molecule, model, conformers)

        for bond, bond_order in zip(molecule.bonds, bond_orders):
            bond.fractional_bond_order = bond_order


class AromaticityProvider(ModelProvider, abc.ABC):
    """The base for classes that provide functionality for applying particular
    aromaticity models to molecules.
    """

    @classmethod
    @abc.abstractmethod
    def _find_aromatic(
        cls, molecule: Molecule, model: str
    ) -> Tuple[AnnotatedList[int], AnnotatedList[int]]:
        """The internal implementation of ``find_aromatic``"""
        raise NotImplementedError()

    @classmethod
    def find_aromatic(
        cls, molecule: Molecule, model: str
    ) -> Tuple[AnnotatedList[int], AnnotatedList[int]]:
        """Returns the indices of any aromatic atoms and bonds in the specified molecule
        by applied a specific aromaticity model.

        Parameters
        ----------
        molecule
            The molecule to apply the aromaticity model to.
        model
            The aromaticity model to apply.
        """

        if model not in cls.supported_models():
            raise UnsupportedFeatureError(cls, model, cls.supported_models())

        return cls._find_aromatic(molecule, model)

    @classmethod
    def assign_aromaticity_model(
        cls,
        molecule: Molecule,
        model: str,
    ):
        """Applies a particular of aromaticity model to the specified molecule, and
        directly assigns the aromaticity flags to the molecule in-place.

        Parameters
        ----------
        molecule
            The molecule to apply the aromaticity model to.
        model
            The aromaticity model to apply.
        """

        aromatic_atoms, aromatic_bonds = cls.find_aromatic(molecule, model)

        for i, atom in enumerate(molecule.atoms):
            atom._is_aromatic = i in aromatic_atoms

        for i, bond in enumerate(molecule.bonds):
            bond._is_aromatic = i in aromatic_bonds


class StereochemistryProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for applying particular
    stereochemistry models to molecules.
    """

    @classmethod
    @abc.abstractmethod
    def find_stereocenters(
        cls, molecule: Molecule
    ) -> Tuple[
        AnnotatedDict[int, AtomStereochemistry], AnnotatedDict[int, BondStereochemistry]
    ]:
        """Returns the indices of stereocenters and bonds in the specified molecule as
        well as their assigned stereo.

        .. todo:: Return type should be literal of supported values and not string.

        Parameters
        ----------
        molecule
            The molecule to check for stereocenters.
        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def _perceive_3d_stereochemistry(
        cls, molecule: Molecule, conformer: unit.Quantity
    ) -> Tuple[
        AnnotatedDict[int, AtomStereochemistry], AnnotatedDict[int, BondStereochemistry]
    ]:
        """The internal implementation of ``perceive_stereochemistry``."""
        raise NotImplementedError()

    @classmethod
    def perceive_3d_stereochemistry(
        cls, molecule: Molecule, conformer: unit.Quantity
    ) -> Tuple[
        AnnotatedDict[int, AtomStereochemistry], AnnotatedDict[int, BondStereochemistry]
    ]:
        """Attempts to perceive the stereochemistry of a molecule based on a 3D
        conformer.

        Parameters
        ----------
        molecule
            The molecule to check for stereocenters.
        conformer
            The 3D conformer of the molecule.
        """

        # TODO: validate the conformer shape and inner type.
        return cls._perceive_3d_stereochemistry(molecule, conformer)


class TautomerProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for enumerating the tautomers
    of molecules.
    """

    @classmethod
    @abc.abstractmethod
    def _enumerate_tautomers(
        cls, molecule: Molecule, max_states: int = 20
    ) -> AnnotatedList[Molecule]:
        """The internal implementation of `enumerate_tautomers`."""
        raise NotImplementedError()

    @classmethod
    def enumerate_tautomers(
        cls, molecule: Molecule, max_states: int = 20
    ) -> AnnotatedList[Molecule]:
        """Enumerate the possible tautomers of a current molecule.

        Parameters
        ----------
        molecule
            The molecule whose state should be enumerated.
        max_states
            The maximum number of tautomers to return.

        Returns
        -------
            A list of tautomers *including* the original tautomer.
        """

        tautomers = cls._enumerate_tautomers(molecule, max_states)

        if len(tautomers) == 0:
            # Handle the case where the input molecule does not have any tautomers.
            return [molecule]

        deduplicate_molecules(tautomers, [molecule])

        return [molecule] + tautomers[: max_states - 1]


class StereoisomerProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for enumerating the stereoisomers
    of molecules.
    """

    @classmethod
    @abc.abstractmethod
    def _enumerate_stereoisomers(
        cls,
        molecule: Molecule,
        undefined_only: bool = False,
        max_isomers: int = 20,
    ) -> AnnotatedList[Molecule]:
        """The internal implementation of `enumerate_stereoisomers`."""
        raise NotImplementedError()

    @classmethod
    def enumerate_stereoisomers(
        cls,
        molecule: Molecule,
        undefined_only: bool = False,
        max_states: int = 20,
    ) -> AnnotatedList[Molecule]:
        """Enumerate the possible stereoisomers of a molecule.

        Parameters
        ----------
        molecule
            The molecule whose state should be enumerated.
        undefined_only
            If all stereocenters and bonds should be enumerate or only those with
            undefined stereochemistry.
        max_states: int optional, default=20
            The maximum number of stereoisomers to return.

        Returns
        --------
            A list of stereoisomers *including* the original stereoisomer if all
            stereochemistry is defined or if the molecule has no stereocenters / bonds.
        """

        stereoisomers = cls._enumerate_stereoisomers(
            molecule, undefined_only, max_states
        )

        if len(stereoisomers) == 0:
            # Handle the case where the input molecule does not have any stereoisomers.
            return [molecule]

        deduplicate_molecules(stereoisomers, [molecule])

        # Check to see if the original molecule had undefined stereochemistry.
        undefined_atom_stereochemistry = any(
            atom_old.stereochemistry is None and atom_new.stereochemistry is not None
            for atom_old, atom_new in zip(molecule.atoms, stereoisomers[0].atoms)
        )
        undefined_bond_stereochemistry = any(
            bond_old.stereochemistry is None and bond_new.stereochemistry is not None
            for bond_old, bond_new in zip(molecule.bonds, stereoisomers[0].bonds)
        )

        # If the input molecule had all of the input stereochemistry defined make sure to
        # return it in addition to the newly found stereoisomers.
        if not undefined_atom_stereochemistry and not undefined_bond_stereochemistry:
            stereoisomers = [molecule] + stereoisomers[: max_states - 1]

        return stereoisomers


class SMARTSProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for matching SMARTS patterns
    against molecules.
    """


class SMILESProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for mapping molecules to and from
    SMILES representations.
    """

    @classmethod
    @abc.abstractmethod
    def _molecule_from_smiles(cls, smiles: str) -> Molecule:
        """The internal implementation of molecule_from_smiles"""
        raise NotImplementedError()

    @classmethod
    def molecule_from_smiles(
        cls, smiles: str, aromaticity_model: str = "mdl"
    ) -> Molecule:

        from chemiwrap.registry import current_registry

        aromaticity_provider = current_registry().aromaticity_providers(
            aromaticity_model
        )

        if aromaticity_provider is None:

            raise MissingFeatureProviderError(
                AromaticityProvider.__name__,
                None,
                f"An aromaticity provider that supports the requested "
                f"{aromaticity_model} model could not be found.",
            )

        stereochemistry_provider = current_registry().stereochemistry_provider()

        if stereochemistry_provider is None:

            raise MissingFeatureProviderError(
                StereoisomerProvider.__name__,
                None,
                "A stereochemistry provider must be available when parsing SMILES "
                "patterns.",
            )

        molecule = cls._molecule_from_smiles(smiles)
        aromaticity_provider.assign_aromaticity_model(molecule, aromaticity_model)

        # Ensure that all of the stereochemistry flags on the molecule are properly
        # set, including the undefined stereochemistry flags.
        atom_stereo, bond_stereo = stereochemistry_provider.find_stereocenters(molecule)

        for i, stereo in atom_stereo.items():
            molecule.atoms[i].stereochemistry = stereo

        return molecule

    @classmethod
    @abc.abstractmethod
    def _molecule_to_smiles(
        cls,
        molecule: Molecule,
        isomeric: bool,
        explicit_hydrogens: bool,
        atom_map: Optional[Dict[int, int]],
    ) -> str:
        """The internal implementation of ``molecule_to_smiles``"""
        raise NotImplementedError()

    @classmethod
    def molecule_to_smiles(
        cls,
        molecule: Molecule,
        isomeric: bool = True,
        explicit_hydrogens: bool = True,
        mapped: bool = False,
    ) -> str:
        """Converts a molecule into a SMILES representation.

        Parameters
        ----------
        molecule
            The molecule to convert into a SMILES representation.
        isomeric
            Whether to record the molecules stereochemistry in the returned
            representation.
        explicit_hydrogens
            Whether to explicitly include all hydrogens in the returned
            representation.
        mapped
            Whether to include atom indices in the returned representation. The atom
            indices can be controlled by supplying an atom map into the molecules'
            properties dictionary. If no mapping is passed all atoms will be mapped in
            order.

        Returns
        -------
            The SMILES representation of the input molecule.
        """

        if mapped and not explicit_hydrogens:

            raise InvalidArgumentError(
                "``explicit_hydrogens`` must be true when creating mapped SMILES."
            )

        # if we only want to map specific atoms check for an atom map
        atom_map = molecule._properties.get("atom_map", None)

        if atom_map is not None:

            map_ids = set(atom_map.values())

            # make sure there are no repeated indices
            if len(map_ids) < len(atom_map):
                atom_map = None

            elif 0 in atom_map.values():

                # we need to increment the map index
                for atom, map in atom_map.items():
                    atom_map[atom] = map + 1

        else:

            atom_map = {i: i + 1 for i in range(molecule.n_atoms)}

        return cls._molecule_to_smiles(
            molecule, isomeric, explicit_hydrogens, None if not mapped else atom_map
        )


class FileIOProvider(FeatureProvider, abc.ABC):
    """The base for classes that provide functionality for reading and writing
    representations of molecules on disk.
    """

    @classmethod
    @abc.abstractmethod
    def supported_export_formats(cls) -> List[str]:
        """The file formats that this provider can export molecules to."""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def supported_import_formats(cls) -> List[str]:
        """The file formats that this provider can import molecules from."""
        raise NotImplementedError()
