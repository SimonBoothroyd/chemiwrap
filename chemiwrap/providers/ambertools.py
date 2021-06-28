from distutils.spawn import find_executable
from typing import Dict, List, Optional, Tuple, Union

from openff.toolkit.topology import Molecule
from simtk import unit
from typing_extensions import Literal, get_args

from chemiwrap.providers import (
    BondOrderProvider,
    ChargeProvider,
    DefaultConformerSettings,
)
from chemiwrap.utilities.provenance import AnnotatedList, AnnotatedQuantity

SQMChargeModel = Literal["gasteiger", "am1bcc", "am1-mulliken"]
SQMBondOrderModel = Literal["am1-wiberg"]


class AmberToolsChargeProvider(ChargeProvider):
    """A wrapper around the charge calculation functionalities provided by Antechamber
    and SQM."""

    @classmethod
    def backend(cls) -> Literal["ambertools"]:
        return "ambertools"

    @classmethod
    def is_available(cls) -> bool:
        return (
            find_executable("sqm") is not None
            and find_executable("antechamber") is not None
        )

    @classmethod
    def supported_models(cls) -> Tuple[SQMChargeModel, ...]:
        return get_args(SQMChargeModel)

    @classmethod
    def default_conformer_settings(
        cls,
    ) -> Dict[SQMChargeModel, Optional[DefaultConformerSettings]]:

        return {
            "am1bcc": DefaultConformerSettings(min_conformers=1, default_conformers=1),
            "am1-mulliken": DefaultConformerSettings(
                min_conformers=1, default_conformers=1
            ),
            "gasteiger": None,
        }

    @classmethod
    def _generate_partial_charges(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Union[List[unit.Quantity], AnnotatedList[unit.Quantity]],
    ) -> AnnotatedQuantity:

        raise NotImplementedError()


class AmberToolsBondOrderProvider(BondOrderProvider):
    """A wrapper around the fraction bond order calculation functionalities provided by
    Antechamber and SQM."""

    @classmethod
    def backend(cls) -> Literal["ambertools"]:
        return "ambertools"

    @classmethod
    def is_available(cls) -> bool:
        return (
            find_executable("sqm") is not None
            and find_executable("antechamber") is not None
        )

    @classmethod
    def supported_models(cls) -> Tuple[SQMChargeModel, ...]:
        return get_args(SQMChargeModel)

    @classmethod
    def default_conformer_settings(
        cls,
    ) -> Dict[str, Optional[DefaultConformerSettings]]:
        return {
            "am1bcc": DefaultConformerSettings(min_conformers=1, default_conformers=1),
            "am1-mulliken": DefaultConformerSettings(
                min_conformers=1, default_conformers=1
            ),
            "gasteiger": None,
        }

    @classmethod
    def _generate_fractional_bond_orders(
        cls,
        molecule: Molecule,
        model: str,
        conformers: Union[List[unit.Quantity], AnnotatedList[unit.Quantity]],
    ) -> AnnotatedList[float]:
        pass
