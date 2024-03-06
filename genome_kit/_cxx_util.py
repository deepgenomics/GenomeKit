# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


def mock(arg):
    """Marks a class, attribute, or method as being a 'mock'.

    When Python-defined class `A` inherits from C++-defined class `B`,
    the attributes and methods `B` are only partially visible to the IDE.

    The `mock` function allows the definition of `A` to mock the attributes
    and methods of `B`, making them accessible the IDE (for autocomplete)
    and to automatic doc extraction (e.g. sphinx).

    For example, the C++-defined :py:class:`genome_kit._cxx.GeneTable` implements
    `__getitem__`, which returns an instance of the Python-defined type
    :py:class:`genome_kit.Gene`. The Python-defined class
    :py:class:`genome_kit.GeneTable` is defined as::

        @_cxx.register
        class GeneTable(_cxx.GeneTable):

            @mock
            def __getitem__(self, index):
                "Returns a Gene object"
                return mock_result(Gene)

    The IDE and autodocs ignore the ``@mock`` decorator, so they treat
    this definition at face value and do static analysis, propagating
    the return type (:py:class:`genome_kit.Gene`) when doing autocomplete::

        for gene in genome.genes:      # genes is of type GeneTable
            gene.<autocomplete>        # gene is of type Gene

    When a class is registered with the C++ backend, all mock
    properties and methods are deleted from the type, so that all
    requests for those attributes will fall through directly to
    the C++ implementation in the base class.

    See the ``PyDeleteMockAttrs`` C++ function for the code that strips
    some mock symbols from each registered type.
    """

    attr = arg
    if isinstance(attr, property):
        attr = attr.fget
    if isinstance(attr, staticmethod):
        attr = attr.__func__
    setattr(attr, "__mock__", True)
    return arg


def mock_unreachable():  # pragma: no cover
    """Prints an unreachable code warning message.

    If a mock attribute or method is ever called, it means there's
    some mismatch between a Python-defined type and a C++-defined type.

    An alternative behaviour would be to raise an exception, but that
    causes headaches when pylint and PyCharm try to warn about
    unreachable code.

    The reason this function is provided separately from `mock_result`
    is that it allows the user to either return nothing::

        def foo(self):
            mock_unreachable()

    or to return an object that provides more complex type hints::

        def foo(self):
            mock_unreachable()
            return [str()]  # type: List[str]

    See also: :py:func:`~genome_kit._cxx_util.mock_result`
    """
    # Make sure we print error message before trying to
    # instantiate the type, because probably instantiation
    # will fail with some error. Don't raise exception
    # or quit(), so we avoid unreachable code warnings
    import os
    if "SPHINXBUILD" in os.environ:
        return
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Mock implementation should be unreachable.")
    print("Name mismatch with C++ attribute/function?")
    print("Call stack:")
    import traceback
    trace = traceback.format_stack()
    print("\n".join(["   " + _ for _ in trace]))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


def mock_result(result_type):  # pragma: no cover
    """Return a mock result of one or more types.

    This function is used to give IDEs a hint during static analysis.

    This function merely calls :py:func:`~genome_kit._cxx.mock_unreachable`
    before returning an instance of `result_type`. It's important to print
    the error message before trying to instantiate the return type.

    See also: :py:func:`~genome_kit._cxx_util.mock`
    """
    mock_unreachable()

    # Do not try to provide more fancy functionality that simply returning the
    # result_type directly. PyCharm gets easily confused by anything more flexible.
    return result_type()


def strip_mock_bases(subtype):
    """Strips the last base classes at runtime from the given type.

    Given a class hierarchy `X->Y->Z` where `X` and `Z` are defined in Python
    but `Y` is defined in C++, then the IDE does not know that `Z` inherits
    from `X`, and so it cannot provide autocomplete for members of `X`.

    This decorator is a hacky way to trick IDEs into statically parsing out
    the inheritance structure, while removing the 'trick' at runtime.

    For example, the class hierarchy we want for `Gene` is simply::

        @_cxx.register
        class Gene(_cxx.Gene):  # Interval -> _cxx.Gene -> Gene
            ...

    However, IDEs will not see that `Gene` inherits from `Interval`.
    So, instead we declare `Gene` as follows::

        @strip_mock_bases
        @_cxx.register
        class Gene(_cxx.Gene, Interval):
            ...

    The IDE will now think that `Gene` inherits `Interval` via multiple
    inheritance, and will IDE autocomplete accordingly, even though at
    runtime we'll get the linear inheritance scheme that we want.
    """
    subtype.__bases__ = subtype.__bases__[:1]
    return subtype


def replace_mock_attr(objtype, obj, name, value):
    """Sets an instance attribute, while also removing mock class
    attribute if it still exists.

    If the `objtype` argument is not specified, it is taken to be `type(obj)`.
    """

    # If the instance's class still has its mock property, delete it.
    if hasattr(objtype, name):
        attr = getattr(objtype, name)
        #if isinstance(attr, property):  # This case not currently used.
        #    attr = attr.fget
        assert hasattr(attr, "__mock__")
        delattr(objtype, name)

    # Set the instance property in its place
    assert not hasattr(obj, name)
    setattr(obj, name, value)
