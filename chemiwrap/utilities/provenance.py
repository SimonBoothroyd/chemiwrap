import copy
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple, TypeVar, Union

from simtk import unit

_T = TypeVar("_T")

_KT = TypeVar("_KT")
_VT = TypeVar("_VT")


class AnnotatedList(Sequence[_T]):
    """An immutable list like collection that has been annotated with provenance."""

    @property
    def provenance(self) -> Dict[str, Any]:
        """The provenance attached to the items in this list."""
        return copy.deepcopy(self.__provenance)

    def __init__(self, iterable: Iterable[_T], provenance: Dict[str, Any]):
        """
        Parameters
        ----------
        iterable
            The items to populate the list with.
        provenance
            The provenance attached to the items.
        """
        self.__list = list(iterable)
        self.__provenance = provenance

    def __contains__(self, item: _T) -> bool:
        return self.__list.__contains__(item)

    def __iter__(self):
        return self.__list.__iter__()

    def __reversed__(self):
        return self.__list.__reversed__()

    def __len__(self) -> int:
        return self.__list.__len__()

    def __getitem__(self, i: int) -> _T:
        return self.__list.__getitem__(i)

    def __repr__(self):
        return self.__list.__repr__()

    def __str__(self):
        return self.__list.__str__()


class AnnotatedDict(Mapping[_KT, _VT]):
    """An immutable dictionary like collection that has been annotated with provenance."""

    @property
    def provenance(self) -> Dict[str, Any]:
        """The provenance attached to the items in this list."""
        return copy.deepcopy(self.__provenance)

    def __init__(
        self,
        iterable: Union[Iterable[Tuple[_KT, _VT]], Mapping[_KT, _VT]],
        provenance: Dict[str, Any],
    ):
        """
        Parameters
        ----------
        iterable
            The items to populate the dictionary with.
        provenance
            The provenance attached to the items.
        """
        self.__dict = dict(iterable)
        self.__provenance = provenance

    def __contains__(self, key: _KT) -> bool:
        return self.__dict.__contains__(key)

    def __iter__(self):
        return self.__dict.__iter__()

    def __reversed__(self):
        return self.__dict.__reversed__()

    def __len__(self) -> int:
        return self.__dict.__len__()

    def __getitem__(self, key: _KT) -> _VT:
        return self.__dict.__getitem__(key)

    def __repr__(self):
        return self.__dict.__repr__()

    def __str__(self):
        return self.__dict.__str__()


class AnnotatedQuantity(unit.Quantity):
    """A unit wrapped quantity that has been annotated with provenance."""

    @property
    def provenance(self) -> Dict[str, Any]:
        """The provenance attached to the items in this list."""
        return copy.deepcopy(self.__provenance)

    def __init__(self, value, unit, provenance: Dict[str, Any]):
        """
        Parameters
        ----------
        provenance
            The provenance attached to the wrapped quantity.
        """
        super(AnnotatedQuantity, self).__init__(value, unit)
        self.__provenance = provenance
