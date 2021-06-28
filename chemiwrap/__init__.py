"""
chemiwrap

A proof of concept refactoring of the OpenFF toolkit registries and wrappers.
"""

from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
