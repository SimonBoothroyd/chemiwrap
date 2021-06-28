from openff.toolkit.topology import Molecule
from simtk import unit

from chemiwrap.registry import FeatureRegistry, current_registry


def main():

    # Using the 'composition' approach we can now easily change the aromaticity
    # model used when parsing SMILES. The aromaticity model will be applied by
    # the provider with the highest priority that supports the model.
    molecule_a = current_registry().smiles_provider().molecule_from_smiles(
        smiles="C1=COC=C1", aromaticity_model="mdl"
    )
    molecule_b = current_registry().smiles_provider().molecule_from_smiles(
        smiles="C1=COC=C1", aromaticity_model="openeye"
    )

    # Same with the library used for stereochemistry perception.
    chiral_amine_smiles = "N(C)(CC)(CCC)"

    with FeatureRegistry(stereochemistry_providers=["openeye"]) as feature_registry:

        oe_stereo = feature_registry.smiles_provider().molecule_from_smiles(
            smiles=chiral_amine_smiles
        )

    with FeatureRegistry(stereochemistry_providers=["rdkit"]) as feature_registry:

        rd_stereo = feature_registry.smiles_provider().molecule_from_smiles(
            smiles=chiral_amine_smiles
        )

    # Create a molecule using the default SMILES provider
    molecule = Molecule.from_smiles("CCCC")

    with FeatureRegistry(
        conformer_providers=["openeye"],
        charge_providers=["openeye"],
    ) as feature_registry:

        # We can access features using the local context manager variable.
        print(feature_registry.providers)

        # Or the 'global' current registry.
        conformers = (
            current_registry()
            .conformer_provider()
            .generate_conformers(
                molecule, max_conformers=10, rms_tolerance=1.0 * unit.angstrom
            )
        )
        print(conformers.provenance)

        # Give RDKit priority for the providers it supports.
        with FeatureRegistry.using_backends("rdkit", "openeye"):

            print(current_registry().providers)

            # Find a charge provider that can provide both AM1 and AM1BCC charges
            # models.
            charge_provider = current_registry().charge_provider(
                ["am1bcc", "am1-mulliken"]
            )

            am1_charges = charge_provider.generate_partial_charges(
                molecule, "am1-mulliken", conformers
            )
            print(am1_charges.provenance["conformers"].provenance)

            am1bcc_charges = charge_provider.generate_partial_charges(
                molecule, "am1bcc"
            )
            print(am1bcc_charges.provenance["conformers"].provenance)


if __name__ == "__main__":
    main()
