"""Code generation."""
from uflx.integrals import AbstractIntegral


def generate(form: AbstractIntegral) -> str:
    """Generate code."""

    print(form)
    print(form.integrand)
    print(form.integrand.arguments)
    print(form.measure)

    code = (
        "void tabulate_tensor_f64(\n"
        "    double* restrict A, const double* restrict w, const double* restrict c,\n"
        "    const double* restrict coordinate_dofs,\n"
        "    const int* restrict entity_local_index,\n"
        "    const uint8_t* restrict quadrature_permutation, void* custom_data\n"
        ") {\n"
        "  static const double weights[3] = {0.1666666666666667, 0.1666666666666667, 0.1666666666666667};\n"
        "  static const double basis_values[1][1][3][3] = {{{{0.6666666666666667, 0.1666666666666666, 0.1666666666666667},\n"
        "    {0.1666666666666667, 0.1666666666666666, 0.6666666666666665},\n"
        "    {0.1666666666666668, 0.6666666666666665, 0.1666666666666667}}}};\n"
        "  static const double FE1_C0_D10_Q48e[1][1][1][3] = {{{{-1.0, 1.0, 0.0}}}};\n"
        "  static const double FE1_C1_D01_Q48e[1][1][1][3] = {{{{-1.0, 0.0, 1.0}}}};\n"
        "  double J0_c0 = 0.0;\n"
        "  double J0_c3 = 0.0;\n"
        "  double J0_c1 = 0.0;\n"
        "  double J0_c2 = 0.0;\n"
        "  for (int ic = 0; ic < 3; ++ic)\n"
        "  {\n"
        "    J0_c0 += coordinate_dofs[ic * 3] * FE1_C0_D10_Q48e[0][0][0][ic];\n"
        "    J0_c3 += coordinate_dofs[ic * 3 + 1] * FE1_C1_D01_Q48e[0][0][0][ic];\n"
        "    J0_c1 += coordinate_dofs[ic * 3] * FE1_C1_D01_Q48e[0][0][0][ic];\n"
        "    J0_c2 += coordinate_dofs[ic * 3 + 1] * FE1_C0_D10_Q48e[0][0][0][ic];\n"
        "  }\n"
        "  double jdet = J0_c0 * J0_c3 - J0_c1 * J0_c2;\n"
        "  double abs_jdet = fabs(jdet);\n"
        "  for (int iq = 0; iq < 3; ++iq)\n"
        "  {\n"
        "    double fw0 = abs_jdet * weights[iq];\n"
        "    double temp_0[3] = {0};\n"
        "    for (int j = 0; j < 3; ++j)\n"
        "      temp_0[j] = fw0 * basis_values[0][0][iq][j];\n"
        "    for (int j = 0; j < 3; ++j)\n"
        "      for (int i = 0; i < 3; ++i)\n"
        "        A[3 * i + j] += basis_values[0][0][iq][i] * temp_0[j];\n"
        "  }\n"
        "}\n"
    )

    signatures = {form: "void tabulate_tensor_f64(double* restrict, const double* restrict, const double* restrict, const double* restrict, const int* restrict, const uint8_t* restrict, void*);"}

    return code, signatures
