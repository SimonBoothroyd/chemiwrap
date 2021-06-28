from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy
from openff.toolkit.topology import Molecule
from simtk import unit

from chemiwrap.utilities.typing import AtomStereochemistry, BondStereochemistry

if TYPE_CHECKING:

    try:
        from openeye.oechem import OEAtomBase, OEBondBase, OEMol
    except ModuleNotFoundError:
        OEAtomBase, OEBondBase, OEMol = Any, Any, Any


def to_openeye(molecule: Molecule) -> "OEMol":
    """Create an OpenEye molecule from an OpenFF molecule.

    Notes
    -----
    * This method no longer assigns aromaticity. Either molecules should not store
      their aromaticity, or the defined values *must* be respected. [DISCUSSION]

    Parameters
    ----------
    molecule
        The molecule to convert to an OEMol

    Returns
    -------
        The OpenEye representation of the molecule.
    """

    from openeye import oechem

    oe_mol = oechem.OEMol()
    oe_mol.SetTitle(molecule.name)

    map_atoms = {}  # {off_idx : oe_idx}

    oe_atoms = list()

    for i, atom in enumerate(molecule.atoms):

        oe_atom = oe_mol.NewAtom(atom.atomic_number)
        oe_atom.SetFormalCharge(
            atom.formal_charge.value_in_unit(unit.elementary_charge)
        )

        oe_atom.SetAromatic(atom.is_aromatic)

        oe_atom.SetData("name", atom.name)
        oe_atom.SetPartialCharge(float("nan"))

        oe_atoms.append(oe_atom)

        map_atoms[i] = oe_atom.GetIdx()

    oe_bonds = list()

    for bond in molecule.bonds:

        atom1_index = bond.atom1_index
        atom2_index = bond.atom2_index

        oe_bond = oe_mol.NewBond(oe_atoms[atom1_index], oe_atoms[atom2_index])
        oe_bond.SetOrder(bond.bond_order)

        oe_bond.SetAromatic(bond.is_aromatic)

        if not (bond.fractional_bond_order is None):
            oe_bond.SetData("fractional_bond_order", bond.fractional_bond_order)

        oe_bonds.append(oe_bond)

    # Set the stereochemistry now that all connectivity is in place
    for atom, oe_atom in zip(molecule.atoms, oe_atoms):

        if not atom.stereochemistry or atom.stereochemistry == "":
            continue

        # Set arbitrary initial stereochemistry
        neighbours = [*oe_atom.GetAtoms()]

        oe_atom.SetStereo(
            neighbours, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right
        )
        oe_cip_stereo = oechem.OEPerceiveCIPStereo(oe_mol, oe_atom)

        if (
            oe_cip_stereo == oechem.OECIPAtomStereo_S and atom.stereochemistry == "S"
        ) or (
            oe_cip_stereo == oechem.OECIPAtomStereo_R and atom.stereochemistry == "R"
        ):
            continue

        oe_atom.SetStereo(
            neighbours, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left
        )

    for bond, oe_bond in zip(molecule.bonds, oe_bonds):

        if not bond.stereochemistry or bond.stereochemistry == "":
            continue

        # Set arbitrary initial stereochemistry
        oe_atom_1 = oe_atoms[bond.atom1_index]
        oe_atom_2 = oe_atoms[bond.atom2_index]

        oe_atom_1_neighbor = [n for n in oe_atom_1.GetAtoms() if n != oe_atom_2][0]
        oe_atom_2_neighbor = [n for n in oe_atom_2.GetAtoms() if n != oe_atom_1][0]

        oe_bond.SetStereo(
            [oe_atom_1_neighbor, oe_atom_2_neighbor],
            oechem.OEBondStereo_CisTrans,
            oechem.OEBondStereo_Cis,
        )
        oe_cip_stereo = oechem.OEPerceiveCIPStereo(oe_mol, oe_bond)

        if (
            oe_cip_stereo == oechem.OECIPBondStereo_E
            and bond.stereochemistry == "E"
            or oe_cip_stereo == oechem.OECIPBondStereo_Z
            and bond.stereochemistry == "Z"
        ):
            continue

        oe_bond.SetStereo(
            [oe_atom_1_neighbor, oe_atom_2_neighbor],
            oechem.OEBondStereo_CisTrans,
            oechem.OEBondStereo_Trans,
        )

    # Retain any present conformations
    if molecule.n_conformers != 0:
        oe_mol.DeleteConfs()

    for conformer in [] if not molecule.conformers else molecule.conformers:

        oe_mol.NewConf(
            oechem.OEFloatArray(conformer.value_in_unit(unit.angstrom).flatten())
        )

    # Retain charges if present. All atoms are initialized above with a partial charge
    # of NaN.
    if molecule.partial_charges is not None:

        oe_indexed_charges = numpy.zeros(shape=molecule.n_atoms, dtype=numpy.float64)

        for off_idx, charge in enumerate(molecule._partial_charges):

            charge_unitless = charge.value_in_units(unit.elementary_charge)
            oe_indexed_charges[map_atoms[off_idx]] = charge_unitless

        for oe_atom in oe_mol.GetAtoms():

            oe_idx = oe_atom.GetIdx()
            oe_atom.SetPartialCharge(oe_indexed_charges[oe_idx])

    # Retain properties if present
    for key, value in molecule.properties.items():
        oechem.OESetSDData(oe_mol, str(key), str(value))

    return oe_mol


def from_openeye(oe_mol: "OEMol") -> Molecule:
    """Create an OpenFF molecule from an OpenEye molecule.

    Notes
    -----
    * This method no longer assigns aromaticity. The correct aromaticity should
      be applied directly to the OE mol and not guessed by this function.
      [DISCUSSION]

    Parameters
    ----------
    oe_mol
        An OpenEye molecule

    Returns
    -------
        An OpenFF molecule
    """

    from openeye import oechem

    # Copy the input molecule so it doesn't accidentally get perturbed.
    oe_mol = oechem.OEMol(oe_mol)

    if oechem.OEHasImplicitHydrogens(oe_mol):
        oechem.OEAddExplicitHydrogens(oe_mol)

    molecule = Molecule()
    molecule.name = oe_mol.GetTitle()

    for data_pair in oechem.OEGetSDDataPairs(oe_mol):
        molecule.properties[data_pair.GetTag()] = data_pair.GetValue()

    map_atoms = dict()  # {oe_mol_idx: molecule_idx}
    atom_mapping = {}

    for oe_atom in oe_mol.GetAtoms():

        oe_idx = oe_atom.GetIdx()
        map_id = oe_atom.GetMapIdx()

        atomic_number = oe_atom.GetAtomicNum()
        formal_charge = oe_atom.GetFormalCharge() * unit.elementary_charge

        is_aromatic = oe_atom.IsAromatic()

        stereochemistry = atom_cip_stereochemistry(oe_mol, oe_atom)
        name = "" if not oe_atom.HasData("name") else oe_atom.GetData("name")

        atom_index = molecule._add_atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            stereochemistry=stereochemistry,
            name=name,
        )

        map_atoms[oe_idx] = atom_index
        atom_mapping[atom_index] = map_id

    # If we have a full / partial atom map add it to the molecule.
    if {*atom_mapping.values()} != {0}:

        molecule.properties["atom_map"] = {
            idx: map_idx for idx, map_idx in atom_mapping.items() if map_idx != 0
        }

    for oe_bond in oe_mol.GetBonds():

        atom1_index = map_atoms[oe_bond.GetBgnIdx()]
        atom2_index = map_atoms[oe_bond.GetEndIdx()]

        bond_order = oe_bond.GetOrder()
        is_aromatic = oe_bond.IsAromatic()

        stereochemistry = bond_cip_stereochemistry(oe_mol, oe_bond)

        fractional_bond_order = (
            None
            if not oe_bond.HasData("fractional_bond_order")
            else oe_bond.GetData("fractional_bond_order")
        )

        molecule._add_bond(
            atom1_index,
            atom2_index,
            bond_order,
            is_aromatic=is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order,
        )

    for conformer in extract_oe_conformers(oe_mol):
        molecule._add_conformer(conformer)

    partial_charges_by_index = {
        map_atoms[oe_atom.GetIdx()]: oe_atom.GetPartialCharge()
        for oe_atom in oe_mol.GetAtoms()
    }
    partial_charges = numpy.array(
        [partial_charges_by_index[i] for i in range(oe_mol.NumAtoms())]
    )

    molecule.partial_charges = (
        None
        if numpy.isnan(partial_charges).all()
        else numpy.array(partial_charges) * unit.elementary_charge
    )

    return molecule


def atom_cip_stereochemistry(
    oe_mol: "OEMol", oe_atom: "OEAtomBase"
) -> Optional[AtomStereochemistry]:
    """Determine CIP stereochemistry (R/S) for the specified atom

    Parameters
    ----------
    oe_mol
        The molecule of interest
    oe_atom
        The atom whose stereochemistry is to be computed

    Returns
    -------
        'R', 'S', or '' if not a chiral atom stereochemistry, or ``None`` is
        unspecified.
    """
    from openeye import oechem

    if oe_atom.IsChiral() and not oe_atom.HasStereoSpecified():
        return None
    elif not oe_atom.IsChiral():
        return ""

    cip = oechem.OEPerceiveCIPStereo(oe_mol, oe_atom)

    if cip == oechem.OECIPAtomStereo_S:
        return "S"
    elif cip == oechem.OECIPAtomStereo_R:
        return "R"
    elif cip == oechem.OECIPAtomStereo_NotStereo:
        return ""

    # [DISCUSSION] - do not leave hanging if, elif statements.
    raise NotImplementedError()


def bond_cip_stereochemistry(
    oe_mol: "OEMol", oe_bond: "OEBondBase"
) -> Optional[BondStereochemistry]:
    """
    Determine CIP stereochemistry (E/Z) for the specified bond

    Parameters
    ----------
    oe_mol
        The molecule of interest
    oe_bond
        The bond whose stereochemistry is to be computed

    Returns
    -------
        'E', 'Z', or '' if not a chiral atom stereochemistry, or ``None`` is
        unspecified.
    """

    from openeye import oechem

    if not oe_bond.HasStereoSpecified():
        return None

    cip = oechem.OEPerceiveCIPStereo(oe_mol, oe_bond)

    if cip == oechem.OECIPBondStereo_E:
        return "E"
    elif cip == oechem.OECIPBondStereo_Z:
        return "Z"
    elif cip == oechem.OECIPBondStereo_NotStereo:
        return ""

    # [DISCUSSION] - do not leave hanging if, elif statements.
    raise NotImplementedError()


def extract_oe_stereochemistry(
    molecule: Molecule, oe_mol: "OEMol"
) -> Tuple[Dict[int, AtomStereochemistry], Dict[int, BondStereochemistry]]:
    """Extracts the CIP stereochemistry of each atom and bond in a OE molecule."""

    atom_stereo = {
        oe_atom.GetIdx(): atom_cip_stereochemistry(oe_mol, oe_atom)
        for oe_atom in oe_mol.GetAtoms()
    }

    bond_stereo_tuples = {
        tuple(
            sorted([oe_bond.GetBgnIdx(), oe_bond.GetEndIdx()])
        ): bond_cip_stereochemistry(oe_mol, oe_bond)
        for oe_bond in oe_mol.GetBonds()
    }
    bond_stereo = {
        i: bond_stereo_tuples[tuple(sorted([bond.atom1_index, bond.atom2_index]))]
        for i, bond in enumerate(molecule.bonds)
    }

    return atom_stereo, bond_stereo


def extract_oe_conformers(oe_mol: "OEMol") -> List[unit.Quantity]:
    """Returns all of the conformers from an RDKit molecule as a list of unit
    wrapped numpy arrays with shape=(n_atoms, 3).
    """

    conformers = []

    for oe_conformer in [] if not hasattr(oe_mol, "GetConfs") else oe_mol.GetConfs():

        conformer = numpy.zeros((oe_mol.NumAtoms(), 3))

        for atom_index, coordinates in oe_conformer.GetCoords().items():
            conformer[atom_index, :] = coordinates

        conformers.append(conformer * unit.angstrom)

    return conformers


def assign_oe_conformers(oe_mol: "OEMol", conformers: List[unit.Quantity]):
    """Assign a set of conformers to an OE molecule, overwriting any existing
    ones."""

    from openeye import oechem

    oe_mol.DeleteConfs()
    assert len(conformers) > 0, "at least one conformer must be provided/"

    for conformer in conformers:

        oe_mol.NewConf(
            oechem.OEFloatArray(conformer.value_in_unit(unit.angstrom).flatten())
        )
