function [apparentResistivity, phase] = modelagem1DMT(modelr, modelt, frequency);

%for i = 1 : length(frequencies)
%    frequency = frequencies(i);
%    [apparentResistivities, phases] = modelMT(resistivities, thicknesses, frequency);
%    apparentResistivity(i) = apparentResistivities;
%    phase(i) = phases;
%end

for i = 1 : length(frequency)
    frequencia = frequency(i);
    [apparentResistivities, phases] = modelMT(modelr, modelt, frequencia);
    apparentResistivity(i) = apparentResistivities;
    phase(i) = phases;
end

%rhoa = [rhoa phase]


