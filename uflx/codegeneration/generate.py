"""Code generation."""
from uflx.integrals import AbstractIntegral


def generate(form: AbstractIntegral) -> str:
    """Generate code."""

    print(form)
    print(form.integrand)
    print(form.integrand.arguments)
    print(form.measure)

    return ""
