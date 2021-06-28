from typing import TYPE_CHECKING, Iterable, Optional, Type

if TYPE_CHECKING:
    from simtk import unit


class MissingFeatureProviderError(ValueError):
    """An exception raised when a feature provider is required but none are
    available."""

    def __init__(
        self,
        feature_provider: str,
        required_features: Optional[Iterable[str]] = None,
        addition_information: Optional[str] = None,
    ):

        message = (
            f"The current feature registry does not provide a {feature_provider} "
            f"provider"
        )

        if required_features is not None:
            message = f"{message} that supports {', '.join(required_features)}."
        else:
            message = f"{message}."

        if addition_information is not None:
            message = f"{message} {addition_information}"

        super(MissingFeatureProviderError, self).__init__(message)


class UnsupportedFeatureError(ValueError):
    """An exception raised when a requested feature is unavailable."""

    def __init__(
        self,
        class_type: Type[object],
        feature: str,
        allowed_features: Optional[Iterable[str]] = None,
    ):

        message = f"{feature} is not supported by the {class_type.__name__} class."

        if allowed_features is not None:
            message = (
                f"{message} The supported values are {', '.join(allowed_features)}"
            )
        else:
            message = f"{message} No specific features are supported by this class."

        super(UnsupportedFeatureError, self).__init__(message)


class InvalidArgumentError(ValueError):
    """An exception raised when an argument provided to a function is invalid."""


class IncompatibleUnitsError(ValueError):
    """An exception raised when an argument is provided with incompatible units."""

    def __init__(
        self,
        argument: str,
        expected_units: "unit.Unit",
        provided_units: "unit.Unit",
    ):
        super(IncompatibleUnitsError, self).__init__(
            f"``{argument}`` has units of {provided_units}, while units of "
            f"{expected_units} were expected."
        )


class ChargeCalculationError(ValueError):
    def __init__(self, addition_information: Optional[str] = None):

        message = "An unhandled error occurred during charge calculation."

        if addition_information is not None:
            message = f"{message} {addition_information}"

        super().__init__(message)


# class ToolkitUnavailableException(MessageException):
#     """The requested toolkit is unavailable."""
#
#     # TODO: Allow toolkit to be specified and used in formatting/printing exception.
#
#
# class LicenseError(ToolkitUnavailableException):
#     """This function requires a license that cannot be found."""
#
#
# class InvalidToolkitError(MessageException):
#     """A non-toolkit object was received when a toolkit object was expected"""
#
#
# class InvalidToolkitRegistryError(MessageException):
#     """An object other than a ToolkitRegistry or toolkit wrapper was received"""
#
#
# class UndefinedStereochemistryError(MessageException):
#     """A molecule was attempted to be loaded with undefined stereochemistry"""
#
#
# class GAFFAtomTypeWarning(RuntimeWarning):
#     """A warning raised if a loaded mol2 file possibly uses GAFF atom types."""
#
#
# class ChargeMethodUnavailableError(ValueError):
#     """A toolkit does not support the requested partial_charge_method combination"""
#
#
# class IncorrectNumConformersError(MessageException):
#     """The requested partial_charge_method expects a different number of conformers than was provided"""
#
#
# class IncorrectNumConformersWarning(Warning):
#     """The requested partial_charge_method expects a different number of conformers than was provided"""
#
#
# class InvalidIUPACNameError(MessageException):
#     """Failed to parse IUPAC name"""
