import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy
from openff.toolkit.topology import Molecule
from simtk import unit
from typing_extensions import Literal, get_args

from chemiwrap.exceptions import ChargeCalculationError
from chemiwrap.providers import (
    AromaticityProvider,
    AtomStereochemistry,
    BondStereochemistry,
    ChargeProvider,
    ConformerProvider,
    DefaultConformerSettings,
    FileIOProvider,
    SMARTSProvider,
    SMILESProvider,
    StereochemistryProvider,
)
from chemiwrap.utilities.openeye import (
    assign_oe_conformers,
    extract_oe_conformers,
    extract_oe_stereochemistry,
    from_openeye,
    to_openeye,
)
from chemiwrap.utilities.provenance import (
    AnnotatedDict,
    AnnotatedList,
    AnnotatedQuantity,
)

_logger = logging.getLogger(__name__)


OEChargeModel = Literal[
    "mmff94", "am1-mulliken", "am1elf10", "am1bcc", "am1bccnosymspt", "am1bccelf10"
]
OEBondOrderModel = Literal["am1-wiberg"]

OEAromaticityModel = Literal["openeye", "mdl"]


class OpenEyeConformerProvider(ConformerProvider):
    """A wrapper around OpenEye's OMEGA conformer generation tool."""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:

        try:
            from openeye import oechem, oeomega
        except ModuleNotFoundError:
            return False

        return oechem.OEChemIsLicensed() and oeomega.OEOmegaIsLicensed()

    @classmethod
    def _generate_conformers(
        cls,
        molecule: Molecule,
        max_conformers: int,
        rms_tolerance: unit.Quantity = 1.0 * unit.angstrom,
    ) -> AnnotatedList[unit.Quantity]:

        import openeye
        from openeye import oeomega

        oe_mol = to_openeye(molecule)

        omega = oeomega.OEOmega()
        omega.SetMaxConfs(max_conformers)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)
        omega.SetEnergyWindow(15.0)

        omega.SetRMSThreshold(rms_tolerance.value_in_unit(unit.angstrom))

        # Don't generate random stereoisomer if not specified
        omega.SetStrictStereo(True)

        status = omega(oe_mol)

        if status is False:

            omega.SetStrictStereo(False)
            new_status = omega(oe_mol)

            if new_status is False:
                raise Exception("OpenEye Omega conformer generation failed")

        conformers = extract_oe_conformers(oe_mol)

        provenance = {
            "backend": cls.backend(),
            "openeye-version": openeye.__version__,
            "max-conformers": max_conformers,
            "rms-tolerance": str(rms_tolerance),
        }

        return AnnotatedList(conformers, provenance)


class OpenEyeChargeProvider(ChargeProvider):
    """A wrapper around OpenEye's QUACPAC partial charge generation tool."""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:

        try:
            from openeye import oechem, oequacpac
        except ModuleNotFoundError:
            return False

        return oechem.OEChemIsLicensed() and oequacpac.OEQuacPacIsLicensed()

    @classmethod
    def supported_models(cls) -> Tuple[OEChargeModel, ...]:
        return get_args(OEChargeModel)

    @classmethod
    def default_conformer_settings(
        cls,
    ) -> Dict[OEChargeModel, Optional[DefaultConformerSettings]]:

        return {
            "mmff94": None,
            "am1-mulliken": {
                "min_conformers": 1,
                "default_conformers": 1,
            },
            "am1elf10": {
                "min_conformers": 1,
                "default_conformers": 500,
            },
            "am1bcc": {
                "min_conformers": 1,
                "default_conformers": 1,
            },
            "am1bccnosymspt": {
                "min_conformers": 1,
                "default_conformers": 1,
            },
            "am1bccelf10": {
                "min_conformers": 1,
                "default_conformers": 500,
            },
        }

    @classmethod
    def _generate_partial_charges(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Union[List[unit.Quantity], AnnotatedList[unit.Quantity]],
    ) -> AnnotatedQuantity:

        import openeye
        from openeye import oechem, oequacpac

        # [DISCUSSION] - the toolkit does not by default symmetrize vanilla
        #                AM1 partial charges.

        oe_charge_method = {
            "mmff94": oequacpac.OEMMFF94Charges(),
            "am1-mulliken": oequacpac.OEAM1Charges(),
            "am1elf10": oequacpac.OEELFCharges(
                oequacpac.OEAM1Charges(optimize=True, symmetrize=True), 10
            ),
            "am1bcc": oequacpac.OEAM1BCCCharges(),
            "am1bccnosymspt": oequacpac.OEAM1Charges(optimize=False, symmetrize=False),
            "am1bccelf10": oequacpac.OEAM1BCCELF10Charges(),
        }[model]

        error_stream = oechem.oeosstream()

        oechem.OEThrow.SetOutputStream(error_stream)
        oechem.OEThrow.Clear()

        oe_mol = to_openeye(molecule)
        assign_oe_conformers(oe_mol, conformers)

        status = oequacpac.OEAssignCharges(oe_mol, oe_charge_method)

        oechem.OEThrow.SetOutputStream(oechem.oeerr)  # restoring to original state

        # This logic handles errors encountered in #34, which can occur when using
        # ELF10 conformer selection
        if (
            not status
            and "SelectElfPop: issue with removing trans COOH conformers"
            in (error_stream.str().decode("UTF-8"))
        ):

            assert model in ["am1elf10", "am1bccelf10"]
            model = "am1-mulliken" if model == "am1elf10" else "am1bcc"

            _logger.warning(
                f"Charge assignment involving ELF10 conformer selection failed "
                f"due to a known bug (OpenFF toolkit issue #346). Downgrading to "
                f"{model} charge assignment for this molecule. More information is "
                f"available at "
                f"https://github.com/openforcefield/openff-toolkit/issues/346"
            )

            return cls._generate_partial_charges(molecule, model, conformers)

        if status is False:
            raise ChargeCalculationError(error_stream.str().decode("UTF-8"))

        # TODO: Make sure atom mapping remains constant
        charges_by_index = {
            oe_atom.GetIdx(): oe_atom.GetPartialCharge()
            for oe_atom in oe_mol.GetAtoms()
        }

        charges = numpy.array([charges_by_index[i] for i in range(molecule.n_atoms)])

        provenance = {
            "backend": cls.backend(),
            "openeye-version": openeye.__version__,
            "model": model,
            # TODO: Get the 'exact' conformers used in the calculation.
            "conformers": conformers,
        }

        return AnnotatedQuantity(
            value=charges, unit=unit.elementary_charge, provenance=provenance
        )


# class OpenEyeBondOrderProvider(BondOrderProvider):
#     """A wrapper around OpenEye's QUACPAC partial charge generation tool that provides
#     functionality for also computing fractional bond orders."""
#
#     @classmethod
#     def backend(cls) -> Literal["openeye"]:
#         return "openeye"
#
#     @classmethod
#     def is_available(cls) -> bool:
#         return False
#
#     @classmethod
#     def supported_models(cls) -> Tuple[OEBondOrderModel, ...]:
#         return get_args(OEBondOrderModel)
#
#     @classmethod
#     def default_conformer_settings(
#         cls,
#     ) -> Dict[str, Optional[DefaultConformerSettings]]:
#         pass


class OpenEyeAromaticityProvider(AromaticityProvider):
    """A wrapper around the aromaticity models provided by OpenEye's cheminformatics
    toolkit."""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:

        try:
            from openeye import oechem
        except ModuleNotFoundError:
            return False

        return oechem.OEChemIsLicensed()

    @classmethod
    def supported_models(cls) -> Tuple[OEAromaticityModel, ...]:
        return get_args(OEAromaticityModel)

    @classmethod
    def _find_aromatic(
        cls, molecule: Molecule, model: str
    ) -> Tuple[AnnotatedList[int], AnnotatedList[int]]:

        import openeye
        from openeye import oechem

        oe_mol = to_openeye(molecule)

        oe_models = {
            "openeye": oechem.OEAroModelOpenEye,
            "mdl": oechem.OEAroModelMDL,
        }

        oechem.OEAssignAromaticFlags(oe_mol, oe_models[model], True)

        aromatic_atom_indices = [
            oe_atom.GetIdx() for oe_atom in oe_mol.GetAtoms() if oe_atom.IsAromatic()
        ]

        aromatic_bond_tuples = [
            tuple(sorted([oe_bond.GetBgnIdx(), oe_bond.GetEndIdx()]))
            for oe_bond in oe_mol.GetBonds()
            if oe_bond.IsAromatic()
        ]
        aromatic_bond_indices = [
            i
            for i, bond in enumerate(molecule.bonds)
            if tuple(sorted([bond.atom1_index, bond.atom2_index]))
            in aromatic_bond_tuples
        ]

        provenance = {
            "backend": cls.backend(),
            "openeye-version": openeye.__version__,
            "model": model,
        }

        return (
            AnnotatedList(aromatic_atom_indices, provenance),
            AnnotatedList(aromatic_bond_indices, provenance),
        )


class OpenEyeStereochemistryProvider(StereochemistryProvider):
    """A wrapper around the stereochemistry models provided by OpenEye's cheminformatics
    toolkit."""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:

        try:
            from openeye import oechem
        except ModuleNotFoundError:
            return False

        return oechem.OEChemIsLicensed()

    @classmethod
    def find_stereocenters(
        cls, molecule: Molecule
    ) -> Tuple[
        AnnotatedDict[int, AtomStereochemistry],
        AnnotatedDict[int, BondStereochemistry],
    ]:

        import openeye
        from openeye import oechem

        oe_mol = to_openeye(molecule)
        assert oechem.OEPerceiveChiral(oe_mol)

        atom_stereo, bond_stereo = extract_oe_stereochemistry(molecule, oe_mol)

        provenance = {
            "backend": cls.backend(),
            "openeye-version": openeye.__version__,
        }

        return (
            AnnotatedDict(atom_stereo, provenance),
            AnnotatedDict(bond_stereo, provenance),
        )

    @classmethod
    def _perceive_3d_stereochemistry(
        cls, molecule: Molecule, conformer: unit.Quantity
    ) -> Tuple[
        AnnotatedDict[int, AtomStereochemistry],
        AnnotatedDict[int, BondStereochemistry],
    ]:

        import openeye
        from openeye import oechem

        oe_mol = to_openeye(molecule)

        oe_mol.DeleteConfs()
        oe_mol.NewConf(
            oechem.OEFloatArray(conformer.value_in_unit(unit.angstrom).flatten())
        )

        assert oechem.OEPerceiveChiral(oe_mol)
        assert oechem.OE3DToInternalStereo(oe_mol)

        atom_stereo, bond_stereo = cls._extract_stereo(molecule, oe_mol)

        provenance = {
            "backend": cls.backend(),
            "openeye-version": openeye.__version__,
        }

        return (
            AnnotatedDict(atom_stereo, provenance),
            AnnotatedDict(bond_stereo, provenance),
        )


# class OpenEyeTautomerProvider(TautomerProvider):
#     """"""
#
#     @classmethod
#     def backend(cls) -> Literal["openeye"]:
#         return "openeye"
#
#     @classmethod
#     def is_available(cls) -> bool:
#         return False
#
#     @classmethod
#     def _enumerate_tautomers(
#         cls, molecule: Molecule, max_states: int = 20
#     ) -> AnnotatedList[Molecule]:
#         pass
#
#
# class OpenEyeStereoisomerProvider(StereoisomerProvider):
#     """"""
#
#     @classmethod
#     def backend(cls) -> Literal["openeye"]:
#         return "openeye"
#
#     @classmethod
#     def is_available(cls) -> bool:
#         return False
#
#     @classmethod
#     def _enumerate_stereoisomers(
#         cls, molecule: Molecule, undefined_only: bool = False, max_isomers: int = 20
#     ) -> AnnotatedList[Molecule]:
#         pass


class OpenEyeSMARTSProvider(SMARTSProvider):
    """"""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:
        return False


class OpenEyeSMILESProvider(SMILESProvider):
    """A wrapper around OpenEye's SMILES parsing / exporting functionality."""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:

        try:
            from openeye import oechem
        except ModuleNotFoundError:
            return False

        return oechem.OEChemIsLicensed()

    @classmethod
    def _molecule_from_smiles(cls, smiles: str) -> Molecule:

        from openeye import oechem

        oe_mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oe_mol, smiles)

        assert oechem.OEAddExplicitHydrogens(oe_mol)

        # Set partial charges to None, since they couldn't have been stored in a SMILES
        for oe_atom in oe_mol.GetAtoms():
            oe_atom.SetPartialCharge(float("nan"))

        return from_openeye(oechem.OEMol(oe_mol))

    @classmethod
    def _molecule_to_smiles(
        cls,
        molecule: Molecule,
        isomeric: bool,
        explicit_hydrogens: bool,
        atom_map: Optional[Dict[int, int]],
    ) -> str:

        from openeye import oechem

        oe_mol = to_openeye(molecule)

        smiles_options = (
            oechem.OESMILESFlag_Canonical
            | oechem.OESMILESFlag_Isotopes
            | oechem.OESMILESFlag_RGroups
        )

        if isomeric:

            smiles_options |= (
                oechem.OESMILESFlag_AtomStereo | oechem.OESMILESFlag_BondStereo
            )

        if explicit_hydrogens:
            smiles_options |= oechem.OESMILESFlag_Hydrogens

        if atom_map is not None:

            for oe_atom in oe_mol.GetAtoms():

                if oe_atom.GetIdx() not in atom_map:
                    continue

                oe_atom.SetMapIdx(atom_map[oe_atom.GetIdx()])

            smiles_options |= oechem.OESMILESFlag_AtomMaps

        smiles = oechem.OECreateSmiString(oe_mol, smiles_options)
        return smiles


class OpenEyeFileIOProvider(FileIOProvider):
    """"""

    @classmethod
    def backend(cls) -> Literal["openeye"]:
        return "openeye"

    @classmethod
    def is_available(cls) -> bool:
        return False

    @classmethod
    def supported_export_formats(cls) -> List[str]:
        return []

    @classmethod
    def supported_import_formats(cls) -> List[str]:
        return []
