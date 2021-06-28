from typing import List, Optional

from openff.toolkit.topology import Molecule


def deduplicate_molecules(
    molecules: List[Molecule], exclude: Optional[List[Molecule]] = None
) -> List[Molecule]:
    """De-duplicates a list of molecules based on their canonical, isomeric, explicit-H,
    non-mapped smiles pattern while *retaining the order* of the original list.

    Parameters
    ----------
    molecules
        The molecules to deduplicate.
    exclude
        An additional set of molecules that should be excluded from the returned list.

    Returns
    --------
        The de-duplicated molecule list.
    """

    exclude = [] if exclude is None else exclude

    unique_smiles = set(molecule.to_smiles(mapped=False) for molecule in exclude)
    unique_molecules = []

    for molecule in molecules:

        smiles = molecule.to_smiles(mapped=False)

        if smiles in unique_smiles:
            continue

        unique_smiles.add(smiles)
        unique_molecules.append(molecule)

    return unique_molecules
