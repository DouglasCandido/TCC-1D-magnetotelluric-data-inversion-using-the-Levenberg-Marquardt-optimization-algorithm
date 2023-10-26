% 1D Inversion of Magnetotelluric Data using Levenberg-Marquardt (LM) Method.

% Inspired by:
% Ekinci Y L and Demirci A 2008 J. of Appl. Sciences 8 (22) : 4070-4078.
% URL: https://scialert.net/abstract/?doi=jas.2008.4070.4078

% Created by:
% Douglas Mateus Soares Cândido da Silva.

close all;
clear all;
clc;

%global data;
global frequency;
global rhoa;
%global modelr;
%global modelt;
global lr;
global lt;

% Load model for forward modeling
% modelr = [200 10 70]; % modelo de resistividade verdadeiro (Ohm-m)
% modelt = [200 400]; % modelo de espessura verdadeiro (m)
% mmodel = [modelr modelt]; % modelo verdadeiro

dadosReais = [
28.5 2.315 0.072 57.19 22.95
38.5 2.254 0.0425 58.19 22.95
52.0 2.229 0.0244 61.39 4.96
70.5 2.188 0.021 59.09 4.46
95.5 2.180 0.0164 59.89 5.96
129.0 2.162 0.0173 51.19 22.95
174.6 2.151 0.0287 46.89 22.95
236.2 2.208 0.0328 42.79 2.46
319.6 2.194 0.019 36.89 1.65
432.5 2.299 0.0270 32.00 22.95
585.1 2.338 0.059 44.00 6.37
791.7 2.420 0.0506 32.00 2.46
1071.1 2.405 0.0825 37.59 22.95
1449.2 2.308 0.123 45.29 4.15
1960.7 2.397 0.092 50.09 22.95];

periodoreal = dadosReais(:,1);
roa = 10.^dadosReais(:,2);
desviopadraoroa = 10.^dadosReais(:,3);
faserealgraus = dadosReais(:,4);
faserealradianos = pi*dadosReais(:,4)/180;
desviopadraofase = dadosReais(:,5);
w = diag(1./desviopadraoroa);

% Inversion

% load rhoa_ruido1.mat;

%logFrequencies = -3:0.2:3;
%frequency = 10.^logFrequencies;
frequency = 1./periodoreal;
%frequency = periodoreal;

rhoa = roa';

% initial model
%r = 10.^[2.7 1.75 2.9 2.5]; % initial resistivity (Ohm-m)
%t = [20000 30000 350000]; % iniial thickness (m)
%m = [r,t];

r = [500 500 500 500]; % initial resistivity (Ohm-m)
t = [20000 20000 20000]; % iniial thickness (m)
m = [r,t];

lr = length(r);
lt = length(t);

kr = 10; % convergence tolerance

j = 1; % initial iteration
iteration(j) = 1;
itermax = 15;

r = m(1:lr);
t = m(1+lr:lr+lt);

[rhoa_cal, phase] = modelagem1DMT(r, t, frequency); % inversão feita utilizando o modelo inicial

%[rhoa_cal, phase] = modelagem1DMT(modelr, modelt, frequency); % inversão feita utilizando o modelo verdadeiro

rms_err = norm(w*rhoa_cal'-w*rhoa')/sqrt(length(rhoa));
rms_error(j) = rms_err;

while(rms_err > kr)

    % find optimum lambda value using golden section search
    [lambda] = gss_lm(m,0.001,10);

    % Jacobian matrix (using input: r,t,rhoa_cal)
    [J] = jacobian(rhoa,r,t,rhoa_cal);

    % Levenberg-Marquardt algorithm
    jac = inv(J'*J+lambda*eye(size(J'*J)));

    dm = jac*J'*[rhoa-rhoa_cal]';

    m = m + dm';
    r = m(1:lr); % resistivity
    t = m(1+lr:lr+lt); % thickness

    [rhoa_cal, phase] = modelagem1DMT(r, t, frequency);

    rms_err = norm(w*rhoa_cal'-w*rhoa')/sqrt(length(rhoa));

    j = j + 1;

    rms_error(j) = rms_err;

    lambdavetor(j) = lambda;

    iteration(j) = iteration(j-1)+1;

     if (iteration(j) > itermax)
         break
     end

     save("m_est.mat", "m")
     save("rms_error.mat", "rms_error")
     save("lambdavetor.mat", "lambdavetor")

end



