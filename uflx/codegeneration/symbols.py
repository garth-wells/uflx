"""Symbols to use in generated code."""


class VariableNamer:
    """Class for naming variables to ensure that temporary variables do not have clashing names."""

    def __init__(self):
        """Initalise."""
        self.i = -1
        self.n = -1
        self.fe_i = -1
        self.qr_i = -1

    def variable(self) -> str:
        """Get a new variable name."""
        chars = ["i", "j", "k"]
        self.i += 1
        if self.i == len(chars):
            self.i = 0
            self.n += 1
        if self.n == -1:
            return chars[self.i]
        else:
            return f"{chars[self.i]}{self.n}"

    def finite_element_table(self) -> str:
        """Get a new finite element table name."""
        self.fe_i += 1
        return f"FE{self.fe_i}"

    def quadrature_table(self) -> str:
        """Get a new quadrature weight table name."""
        self.qr_i += 1
        return f"QW{self.qr_i}"


global_variable_namer = VariableNamer()

local_tensor = "A"
coordinate_dofs = "coordinate_dofs"
constants = "c"
coefficients = "w"
entity_local_index = "entity_local_index"
quadrature_permutation = "quadrature_permutation"
custom_data = "custom_data"
