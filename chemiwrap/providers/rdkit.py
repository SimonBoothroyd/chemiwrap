import importlib
from typing import Tuple

from openff.toolkit.topology import Molecule
from simtk import unit
from typing_extensions import Literal

from chemiwrap.providers import ConformerProvider, StereochemistryProvider
from chemiwrap.utilities.provenance import AnnotatedDict, AnnotatedList
from chemiwrap.utilities.rdkit import (
    assign_rd_conformers,
    extract_rd_conformers,
    extract_rd_stereochemistry,
    to_rdkit,
)
from chemiwrap.utilities.typing import AtomStereochemistry, BondStereochemistry

RDChargeModel = Literal["mmff94"]


def _rdkit_available() -> bool:

    try:
        importlib.import_module("rdkit")
    except ImportError:
        return False

    return True


class RDKitConformerProvider(ConformerProvider):
    """A wrapper around RDKit's conformer generation functionality."""

    @classmethod
    def backend(cls) -> Literal["rdkit"]:
        return "rdkit"

    @classmethod
    def is_available(cls) -> bool:
        return _rdkit_available()

    @classmethod
    def _generate_conformers(
        cls,
        molecule: Molecule,
        max_conformers: int,
        rms_tolerance: unit.Quantity = 1.0 * unit.angstrom,
    ) -> AnnotatedList[unit.Quantity]:

        import rdkit
        from rdkit.Chem import AllChem

        rd_mol = to_rdkit(molecule)

        # TODO: This generates way more conformations than omega, given the same
        #       nConfs and RMS threshold. Is there some way to set an energy cutoff
        #       as well?
        AllChem.EmbedMultipleConfs(
            rd_mol,
            numConfs=max_conformers,
            pruneRmsThresh=rms_tolerance.value_in_unit(unit.angstrom),
            randomSeed=1,
        )

        conformers = extract_rd_conformers(rd_mol)

        provenance = {
            "backend": cls.backend(),
            "rdkit-version": rdkit.__version__,
            "max-conformers": max_conformers,
            "rms-tolerance": str(rms_tolerance),
        }

        return AnnotatedList(conformers, provenance)


class RDKitStereochemistryProvider(StereochemistryProvider):
    """A wrapper around the stereochemistry models provided by OpenEye's cheminformatics
    toolkit."""

    @classmethod
    def backend(cls) -> Literal["rdkit"]:
        return "rdkit"

    @classmethod
    def is_available(cls) -> bool:
        return _rdkit_available()

    @classmethod
    def find_stereocenters(
        cls, molecule: Molecule
    ) -> Tuple[
        AnnotatedDict[int, AtomStereochemistry],
        AnnotatedDict[int, BondStereochemistry],
    ]:

        import rdkit
        from rdkit import Chem

        rd_mol = to_rdkit(molecule)
        Chem.AssignStereochemistryFrom3D(rd_mol)

        atom_stereo, bond_stereo = extract_rd_stereochemistry(molecule, rd_mol)

        provenance = {
            "backend": cls.backend(),
            "rdkit-version": rdkit.__version__,
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

        import rdkit
        from rdkit import Chem

        rd_mol = to_rdkit(molecule)
        assign_rd_conformers(rd_mol, [conformer])

        Chem.AssignStereochemistryFrom3D(rd_mol)

        atom_stereo, bond_stereo = extract_rd_stereochemistry(molecule, rd_mol)

        provenance = {
            "backend": cls.backend(),
            "rdkit-version": rdkit.__version__,
        }

        return (
            AnnotatedDict(atom_stereo, provenance),
            AnnotatedDict(bond_stereo, provenance),
        )
