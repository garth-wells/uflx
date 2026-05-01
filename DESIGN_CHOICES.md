# Design choices

Documentation of design decisions that have been made

## Naming
File names are plural to allow (eg) `uflx.domains` to be the module and `uflx.domain` to be
the function defined in domains.py that is imported into __init__.py.

Whenever a function to initialise a class is defined, it's name is the same as the class name
but lowercase, with `_x` wherever the class name has (eg) a capital `X` mid name.

## Default values
Class `__init__` functions have no default values. If default values are required, a function
to initialise the class with defaults should be defined.
