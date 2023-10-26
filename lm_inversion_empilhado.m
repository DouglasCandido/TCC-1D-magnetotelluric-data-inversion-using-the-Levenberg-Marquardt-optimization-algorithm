% 1D Inversion of Magnetotelluric Data using Levenberg-Marquardt (LM) Method.

% Inspired by:
% Ekinci Y L and Demirci A 2008 J. of Appl. Sciences 8 (22) : 4070-4078.
% URL: https://scialert.net/abstract/?doi=jas.2008.4070.4078

% Created by:
% Douglas Mateus Soares Cândido da Silva.

close all;
clear;
clc;

%global data;
global frequency;

global rhoa;

global modelr;
global modelt;

global lr;
global lt;

% Load model for forward modeling
modelr = [300 2500 0.8 3000 2500]; % modelo de resistividade verdadeiro (Ohm-m)
modelt = [200 400 40 500]; % modelo de espessura verdadeiro (m)
mmodel = [modelr modelt]; % modelo verdadeiro

% Inversion
load rhoa_ruido.mat;
load phase_ruido.mat;

logFrequencies = -4:0.1:4;
frequency = 10.^logFrequencies;

rhoa = [rhoa_ruido phase_ruido];

% initial model

r = [500,1500,100,2000,1500]; % initial resistivity (Ohm-m)
t = [200,200,20,200]; % iniial thickness (m)
m = [r,t];

lr = length(r);
lt = length(t);

kr = 10e-2; % convergence tolerance

j = 1; % initial iteration
iteration(j) = 1;
itermax = 50;

r = m(1:lr);
t = m(1+lr:lr+lt);

%[rhoa_cal, phase] = modelagem1DMT(r, t, frequency); % inversão feita utilizando o modelo inicial

%[rhoa_cal, phase] = modelagem1DMT(modelr, modelt, frequency); % inversão feita utilizando o modelo verdadeiro

d_emp_cal = modelagem1DMTEmpilhado(modelr, modelt, frequency);

rms_err = norm(d_emp_cal-rhoa)/sqrt(length(rhoa));
rms_error(j) = rms_err;

%if(rms_err < kr)
%  plot_iterasi(m,rhoa_cal,j);
%end

while(rms_err > kr)

    % find optimum lambda value using golden section search
    [lambda] = gss_lmEmp(m,0.001,10);

    %disp('lambda'); disp(lambda);

    %lamda_plot(j) = lambda;

    % Jacobian matrix (using input: r,t,rhoa_cal)
    [J] = jacobianEmp(rhoa,r,t,d_emp_cal);

    % Levenberg-Marquardt algorithm
    jac = inv(J'*J+lambda*eye(size(J'*J)));

    dm = jac*J'*[rhoa-d_emp_cal]';

    m = m + dm';
    r = m(1:lr); % resistivity
    t = m(1+lr:lr+lt); % thickness

    rhoa_cal = modelagem1DMTEmpilhado(r, t, frequency);

    rms_err = norm(rhoa_cal-rhoa)/sqrt(length(rhoa));

    j = j + 1;

    rms_error(j) = rms_err;

    iteration(j) = iteration(j-1)+1;

     if (iteration(j) > itermax)
         break
     end

     %plot_iterasi(m,rhoa_cal,j);

end
