Data used to implement the TADA_RD method (Moon et al., 2025):

<data required to run the TADA_RD method>

  - genetable_info: Gene-level information required to run the TADA_RD method.
  - ClassDn_rusboost: A classifier trained using the SPARK family-based dataset with the RUSBoost algorithm.
  - ClassDn_underbagging: A classifier trained using the SPARK family-based dataset with the UnderBagging algorithm.
  - accuracy_rusboost: Sensitivity and specificity of ClassDn_rusboost at various threshold levels.
  - accuracy_underbagging: Sensitivity and specificity of ClassDn_underbagging at various threshold levels.

<example datasets to run the TADA_RD method>
  - fu_supplementary_table: An example family-based dataset from the supplementary table 5 of Fu et al. (2022).
  - PROBANDS-simulated-ptv-variants: Simulated case data generated using the method described in Supplementary Material B of Moon et al. (2025).
  - SIBLINGS-simulated-ptv-variants: Simulated control data generated using the method described in Supplementary Material B of Moon et al. (2025).

