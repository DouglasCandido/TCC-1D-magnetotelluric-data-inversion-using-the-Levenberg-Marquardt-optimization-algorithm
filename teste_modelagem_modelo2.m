close all;

clear;

clc;

% Load model for forward modeling
%modelr = [30000 5000 1000]; % modelo de resistividade verdadeiro (Ohm-m)
%modelt = [15 18]; % modelo de espessura verdadeiro (Km)
%modelr = [250 25 100 10 1000]; % modelo de resistividade verdadeiro (Ohm-m)
%modelt = [0.6 0.4 2 0.25]; % modelo de espessura verdadeiro (Km)
modelr = [1000 100 1000]; % modelo de resistividade verdadeiro (Ohm-m)
modelt = [1000 1000]; % modelo de espessura verdadeiro (m)
mmodel = [modelr modelt]; % modelo verdadeiro

% initial model
%r = [10000,1000,100]; % initial resistivity (Ohm-m)
%t = [5,8]; % iniial thickness (Km)
r = [50,50,50]; % initial resistivity (Ohm-m)
t = [300,300]; % iniial thickness (m)
%r = [100, 10, 100]; % initial resistivity (Ohm-m)
%t = [0.1 0.1]; % iniial thickness (Km)
m = [r,t];

%%
logFrequencies = -3:0.2:3;
frequency = 10.^logFrequencies;

% m_est = [52.3100    52.3094    52.3080    52.3061   984.6082     1.0157     1.0155     1.0150     1.0144];
%r_est = [52.3100 52.3094    52.3080    52.3061   984.6082]
%[rhoa_est, phase_est] = modelagem1DMT(r_est, t, frequency)

[rhoa, phase] = modelagem1DMT(modelr, modelt, frequency) % inversão feita utilizando o modelo verdadeiro

rhoa_ruido2 = rhoa + 0.05 * randn(size(rhoa)).*rhoa; % acrescentando o ruído na inversão

phase_ruido2 = phase + 0.05 * randn(size(phase)).*phase; % acrescentando o ruído na fase

[rhoa2, phase2] = modelagem1DMT(r, t, frequency)

save("rhoa_ruido2.mat", "rhoa_ruido2")
save("phase_ruido2.mat", "phase_ruido2")

subplot(1,2,1)
%loglog(frequency,rhoa)
loglog(frequency,rhoa2)
%loglog(frequency,rhoa_est)
hold on
loglog(frequency,rhoa_ruido2,"*")
hold off
subplot(1,2,2)
loglog(frequency,phase)
hold on
loglog(frequency,phase_ruido2,"*")

%[J] = jacobian(rhoa_ruido,modelr,modelt,rhoa);
