% 1D Inversion of Magnetotelluric Data using Levenberg-Marquardt (LM) Method.

% Inspired by:
% Ekinci Y L and Demirci A 2008 J. of Appl. Sciences 8 (22) : 4070-4078.
% URL: https://scialert.net/abstract/?doi=jas.2008.4070.4078

% Created by:
% Douglas Mateus Soares Cândido da Silva.

close all;
clear;
clc;

global frequency;
global rhoa;
global modelr;
global modelt;
global lr;
global lt;

% Load model for forward modeling
modelr = [1000 100 1000]; % modelo de resistividade verdadeiro (Ohm-m)
modelt = [1000 1000]; % modelo de espessura verdadeiro (m)
mmodel = [modelr modelt]; % modelo verdadeiro

% Inversion
load rhoa_ruido2.mat;

logFrequencies = -3:0.2:3;
frequency = 10.^logFrequencies;

rhoa = rhoa_ruido2;

% initial model
r = [50,50,50]; % initial resistivity (Ohm-m)
t = [300,300]; % iniial thickness (m)
m = [r,t];

lr = length(r);
lt = length(t);

kr = 2; % convergence tolerance

j = 1; % initial iteration
iteration(j) = 1;
itermax = 10;

r = m(1:lr);
t = m(1+lr:lr+lt);

[rhoa_cal, phase] = modelagem1DMT(r, t, frequency); % inversão feita utilizando o modelo inicial

rms_err = norm(rhoa_cal-rhoa)/sqrt(length(rhoa));
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

    rms_err = norm(rhoa_cal-rhoa)/sqrt(length(rhoa));

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



