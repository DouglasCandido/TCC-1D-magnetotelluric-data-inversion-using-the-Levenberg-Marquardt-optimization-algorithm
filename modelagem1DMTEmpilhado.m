function dadoEmpilhado = modelagem1DMTEmpilhado(modelr, modelt, frequency)

  [apparentResistivity, phase] = modelagem1DMT(modelr, modelt, frequency);

  dadoEmpilhado = [apparentResistivity phase];


