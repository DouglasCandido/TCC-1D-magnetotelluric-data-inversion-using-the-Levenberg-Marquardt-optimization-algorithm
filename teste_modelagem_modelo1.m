close all;

clear;

clc;

% Load model for forward modeling
modelr = [200 10 70]; % modelo de resistividade verdadeiro (Ohm-m)
modelt = [200 400]; % modelo de espessura verdadeiro (m)
mmodel = [modelr modelt]; % modelo verdadeiro

% initial model
r = [50,50,50]; % initial resistivity (Ohm-m)
t = [300,300]; % iniial thickness (Km)
m = [r,t];

%%

logFrequencies = -3:0.2:3;
frequency = 10.^logFrequencies;

%[rhoa, phase] = modelagem1dMT(r, t, frequency)

[rhoa, phase] = modelagem1DMT(modelr, modelt, frequency) % inversão feita utilizando o modelo verdadeiro

rhoa_ruido1 = rhoa + 0.025 * randn(size(rhoa)).*rhoa; % acrescentando o ruído na inversão

phase_ruido1 = phase + 0.025 * randn(size(phase)).*phase; % acrescentando o ruído na fase

save("rhoa_ruido1.mat", "rhoa_ruido1")
save("phase_ruido1.mat", "phase_ruido1")

subplot(1,2,1)
loglog(frequency,rhoa)
hold on
loglog(frequency,rhoa_ruido1,"*")
hold off
subplot(1,2,2)
loglog(frequency,phase)
hold on
loglog(frequency,phase_ruido1,"*")

%[J] = jacobian(rhoa_ruido,modelr,modelt,rhoa);
