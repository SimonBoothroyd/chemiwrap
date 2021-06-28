import itertools
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy
from openff.toolkit.topology import Molecule
from simtk import unit

from chemiwrap.utilities.typing import AtomStereochemistry, BondStereochemistry

if TYPE_CHECKING:

    try:
        from rdkit.Chem import Atom as RDAtom
        from rdkit.Chem import Bond as RDBond
        from rdkit.Chem import Mol as RDMol
    except ModuleNotFoundError:
        RDAtom, RDBond, RDMol = Any, Any, Any


def to_rdkit(molecule: Molecule) -> "RDMol":
    """Create an RDKit molecule from an OpenFF molecule.

    Notes
    -----
    * This method no longer assigns aromaticity. Either molecules should not store
      their aromaticity, or the defined values *must* be respected. [DISCUSSION]

    Parameters
    ----------
    molecule
        The molecule to convert to an RDMol

    Returns
    -------
        The RDKit representation of the molecule.
    """
    from rdkit import Chem, Geometry

    rd_mol = Chem.RWMol()

    if not (molecule.name is None):
        rd_mol.SetProp("_Name", molecule.name)

    for index, atom in enumerate(molecule.atoms):

        rd_atom = Chem.Atom(atom.atomic_number)
        rd_atom.SetFormalCharge(
            atom.formal_charge.value_in_unit(unit.elementary_charge)
        )
        rd_atom.SetIsAromatic(atom.is_aromatic)
        rd_atom.SetProp("_Name", atom.name)

        rd_index = rd_mol.AddAtom(rd_atom)
        assert index == rd_index

    _bond_order_to_type = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.QUADRUPLE,
        5: Chem.BondType.QUINTUPLE,
        6: Chem.BondType.HEXTUPLE,
        7: Chem.BondType.ONEANDAHALF,
    }

    for bond in molecule.bonds:

        atom_indices = (
            bond.atom1.molecule_atom_index,
            bond.atom2.molecule_atom_index,
        )

        rd_mol.AddBond(*atom_indices)
        rd_bond = rd_mol.GetBondBetweenAtoms(*atom_indices)

        if not (bond.fractional_bond_order is None):
            rd_bond.SetDoubleProp("fractional_bond_order", bond.fractional_bond_order)

        if bond.is_aromatic:
            rd_bond.SetBondType(_bond_order_to_type[1.5])
            rd_bond.SetIsAromatic(True)
        else:
            rd_bond.SetBondType(_bond_order_to_type[bond.bond_order])
            rd_bond.SetIsAromatic(False)

    # TODO: is aromaticity presevered by the sanitization?
    Chem.SanitizeMol(
        rd_mol,
        Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY,
    )

    # No longer forcibly set the _CIPCode property of these atoms... it seems super
    # fragile and in the current implementation of to rdkit -> pull properties from
    # rdmol rather than to rdkit -> calc props -> from rdkit should not cause issues...
    for index, atom in enumerate(molecule.atoms):

        if not atom.stereochemistry or atom.stereochemistry == "":
            continue

        rd_atom = rd_mol.GetAtomWithIdx(index)

        # Set an arbitrary initial stereochemistry and then flip if the corresponding
        # CIP does not match.
        rd_atom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(rd_mol, force=True, cleanIt=True)

        cip = atom_cip_stereochemistry(rd_atom)

        if cip == atom.stereochemistry:
            continue

        rd_atom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
        Chem.AssignStereochemistry(rd_mol, force=True, cleanIt=True)

    # Copy bond stereo info from molecule to rdmol. This is a bit complex in RDKit so
    # for now let's just enumerate stereoisomers until we get one that matches...
    rd_neighbour_bonds = {}

    for bond in molecule.bonds:

        # No need to do anything with bonds without stereochemistry.
        if not bond.stereochemistry or bond.stereochemistry == "":
            continue

        atom_1, atom_2 = bond.atom1, bond.atom2

        atom_1_index = atom_1.molecule_atom_index
        atom_2_index = atom_2.molecule_atom_index

        neighbour_1 = [n for n in atom_1.bonded_atoms if n != atom_2][
            0
        ].molecule_atom_index
        neighbour_2 = [n for n in atom_2.bonded_atoms if n != atom_1][
            0
        ].molecule_atom_index

        rd_neighbour_bonds[
            tuple(sorted([neighbour_1, atom_1_index]))
        ] = rd_mol.GetBondBetweenAtoms(atom_1_index, neighbour_1)
        rd_neighbour_bonds[
            tuple(sorted([neighbour_2, atom_2_index]))
        ] = rd_mol.GetBondBetweenAtoms(atom_2_index, neighbour_2)

    rd_neighbour_bonds = [*rd_neighbour_bonds.values()]
    stereo_matches = False

    for directions in itertools.product(
        [Chem.rdchem.BondDir.ENDDOWNRIGHT, Chem.rdchem.BondDir.ENDUPRIGHT],
        repeat=len(rd_neighbour_bonds),
    ):

        for rd_bond, direction in zip(rd_neighbour_bonds, directions):
            rd_bond.SetBondDir(direction)

        Chem.AssignStereochemistry(rd_mol, cleanIt=True, force=True)

        stereo_matches = True

        for bond in molecule.bonds:

            # No need to do anything with bonds without stereochemistry.
            if not bond.stereochemistry or bond.stereochemistry == "":
                continue

            rd_bond = rd_mol.GetBondBetweenAtoms(bond.atom1_index, bond.atom2_index)
            cip = bond_cip_stereochemistry(rd_bond)

            if cip == bond.stereochemistry:
                continue

            stereo_matches = False
            break

        if stereo_matches:
            break

    assert stereo_matches

    # Set coordinates if we have them
    for conformer in [] if molecule.conformers is None else molecule._conformers:

        rd_conformer = Chem.Conformer()

        for atom_idx in range(molecule.n_atoms):
            x, y, z = conformer[atom_idx, :].value_in_unit(unit.angstrom)
            rd_conformer.SetAtomPosition(atom_idx, Geometry.Point3D(x, y, z))

        rd_mol.AddConformer(rd_conformer, assignId=True)

    # Retain charges, if present
    if not (molecule._partial_charges is None):

        rdk_indexed_charges = numpy.zeros(shape=molecule.n_atoms, dtype=float)

        for atom_idx, charge in enumerate(molecule._partial_charges):
            charge_unitless = charge.value_in_unit(unit.elementary_charge)
            rdk_indexed_charges[atom_idx] = charge_unitless
        for atom_idx, rdk_atom in enumerate(rd_mol.GetAtoms()):
            rdk_atom.SetDoubleProp("PartialCharge", rdk_indexed_charges[atom_idx])

        Chem.CreateAtomDoublePropertyList(rd_mol, "PartialCharge")

    # Retain properties if present
    for name, value in molecule.properties.items():

        if type(value) == str:
            rd_mol.SetProp(name, value)
        elif type(value) == int:
            rd_mol.SetIntProp(name, value)
        elif type(value) == float:
            rd_mol.SetDoubleProp(name, value)
        elif type(value) == bool:
            rd_mol.SetBoolProp(name, value)
        else:
            rd_mol.SetProp(name, str(value))

    rd_mol.UpdatePropertyCache(strict=False)

    return Chem.Mol(rd_mol)


def from_rdkit(rd_mol: "RDMol") -> Molecule:
    """Create an OpenFF molecule from an OpenEye molecule.

    Notes
    -----
    * This method no longer assigns aromaticity. The correct aromaticity should
      be applied directly to the OE mol and not guessed by this function.
      [DISCUSSION]

    Parameters
    ----------
    rd_mol
        An RDKit molecule

    Returns
    -------
        An OpenFF molecule
    """
    from rdkit import Chem

    # Copy the input molecule so it doesn't accidentally get perturbed.
    rd_mol = Chem.AddHs(Chem.Mol(rd_mol), addCoords=True)

    Chem.SanitizeMol(rd_mol, Chem.SANITIZE_ALL)
    Chem.AssignStereochemistry(rd_mol, cleanIt=True, force=True)

    molecule = Molecule()
    molecule.name = ""

    if rd_mol.HasProp("_Name"):
        molecule.name = rd_mol.GetProp("_Name")

    molecule._properties = rd_mol.GetPropsAsDict()

    map_atoms = {}  # {rd_mol_idx: molecule_idx}
    atom_mapping = {}

    for rd_atom in rd_mol.GetAtoms():

        rd_idx = rd_atom.GetIdx()

        try:
            map_id = int(rd_atom.GetProp("_map_idx"))
        except KeyError:
            map_id = rd_atom.GetAtomMapNum()

        atomic_number = rd_atom.GetAtomicNum()
        formal_charge = rd_atom.GetFormalCharge() * unit.elementary_charge

        is_aromatic = rd_atom.GetIsAromatic()

        name = "" if not rd_atom.HasProp("_Name") else rd_atom.GetProp("_Name")

        # If chiral, store the chirality to be set later
        stereochemistry = atom_cip_stereochemistry(rd_atom)

        atom_index = molecule._add_atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            name=name,
            stereochemistry=stereochemistry,
        )
        map_atoms[rd_idx] = atom_index
        atom_mapping[atom_index] = map_id

    # If we have a full / partial atom map add it to the molecule.
    if {*atom_mapping.values()} != {0}:

        molecule.properties["atom_map"] = {
            idx: map_idx for idx, map_idx in atom_mapping.items() if map_idx != 0
        }

    for rd_bond in rd_mol.GetBonds():

        atom1_index = rd_bond.GetBeginAtomIdx()
        atom2_index = rd_bond.GetEndAtomIdx()

        bond_order = int(rd_bond.GetBondTypeAsDouble())
        is_aromatic = rd_bond.GetIsAromatic()

        stereochemistry = bond_cip_stereochemistry(rd_bond)

        fractional_bond_order = (
            None
            if not rd_bond.HasProp("fractional_bond_order")
            else rd_bond.GetDoubleProp("fractional_bond_order")
        )

        molecule._add_bond(
            atom1_index,
            atom2_index,
            bond_order,
            is_aromatic=is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order,
        )

    for conformer in extract_rd_conformers(rd_mol):
        molecule._add_conformer(conformer)

    partial_charges_by_index = {
        map_atoms[rd_atom.GetIdx()]: numpy.nan
        if not rd_atom.HasProp("PartialCharge")
        else rd_atom.GetDoubleProp("PartialCharge")
        for rd_atom in rd_mol.GetAtoms()
    }
    partial_charges = numpy.array(
        [partial_charges_by_index[i] for i in range(rd_mol.GetNumAtoms())]
    )

    molecule.partial_charges = (
        None
        if numpy.isnan(partial_charges).all()
        else numpy.array(partial_charges) * unit.elementary_charge
    )

    return molecule


def atom_cip_stereochemistry(rd_atom: "RDAtom") -> Optional[AtomStereochemistry]:
    """Determine CIP stereochemistry (R/S) for the specified atom

    Parameters
    ----------
    rd_atom
        The atom whose stereochemistry is to be evaluated

    Returns
    -------
        'R', 'S', or '' if not a chiral atom stereochemistry, or ``None`` is
        unspecified.
    """

    if not rd_atom.HasProp("_CIPCode"):
        return ""

    cip = rd_atom.GetProp("_CIPCode")

    if cip not in {"R", "S"}:
        raise NotImplementedError()

    return cip


def bond_cip_stereochemistry(rd_bond: "RDBond") -> Optional[BondStereochemistry]:
    """
    Determine CIP stereochemistry (E/Z) for the specified bond

    Parameters
    ----------
    rd_bond
        The bond whose stereochemistry is to be evaluated

    Returns
    -------
        'E', 'Z', or '' if not a chiral atom stereochemistry, or ``None`` is
        unspecified.
    """

    from rdkit import Chem

    cip = rd_bond.GetStereo()

    if cip == Chem.rdchem.BondStereo.STEREOANY:
        return None

    if cip == Chem.rdchem.BondStereo.STEREOE:
        return "E"
    elif cip == Chem.rdchem.BondStereo.STEREOZ:
        return "Z"
    elif cip == Chem.rdchem.BondStereo.STEREONONE:
        return ""

    raise NotImplementedError()


def extract_rd_stereochemistry(
    molecule: Molecule, rd_mol: "RDMol"
) -> Tuple[Dict[int, AtomStereochemistry], Dict[int, BondStereochemistry]]:
    """Extracts the CIP stereochemistry of each atom and bond in a OE molecule."""

    atom_stereo = {
        rd_atom.GetIdx(): atom_cip_stereochemistry(rd_atom)
        for rd_atom in rd_mol.GetAtoms()
    }

    bond_stereo_tuples = {
        tuple(
            sorted([rd_bond.GetBeginAtomIdx(), rd_bond.GetEndAtomIdx()])
        ): bond_cip_stereochemistry(rd_bond)
        for rd_bond in rd_mol.GetBonds()
    }
    bond_stereo = {
        i: bond_stereo_tuples[tuple(sorted([bond.atom1_index, bond.atom2_index]))]
        for i, bond in enumerate(molecule.bonds)
    }

    return atom_stereo, bond_stereo


def extract_rd_conformers(rd_mol: "RDMol") -> List[unit.Quantity]:
    """Returns all of the conformers from an RDKit molecule as a list of unit
    wrapped numpy arrays with shape=(n_atoms, 3).
    """

    return [
        rd_conformer.GetPositions() * unit.angstrom
        for rd_conformer in rd_mol.GetConformers()
    ]


def assign_rd_conformers(rd_mol: "RDMol", conformers: List[unit.Quantity]):
    """Assign a set of conformers to an RDKit molecule, overwriting any existing
    ones."""

    from rdkit import Chem, Geometry

    rd_mol.RemoveAllConformers()

    for conformer in conformers:

        rd_conformer = Chem.Conformer()

        for atom_idx in range(rd_mol.GetNumAtoms()):

            x, y, z = conformer[atom_idx, :].value_in_unit(unit.angstrom)
            rd_conformer.SetAtomPosition(atom_idx, Geometry.Point3D(x, y, z))

        rd_mol.AddConformer(rd_conformer, assignId=True)
