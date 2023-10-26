close all;
clc;
clear;

%load rhoa_ruido1.mat;
%load phase_ruido1.mat;
load rms_error.mat;

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
desviopadraoroa = dadosReais(:,3);
faserealgraus = dadosReais(:,4);
faserealradianos = pi*dadosReais(:,4)/180;
desviopadraofase = dadosReais(:,5);

rhoa = roa';

frequency = 1./periodoreal;

% Load model for forward modeling
%modelr = [150,150,150,150]; % modelo de resistividade verdadeiro (Ohm-m)
%modelt = [1000,1000,1000]; % modelo de espessura verdadeiro (m)
%mmodel = [modelr modelt]; % modelo verdadeiro

%modelr = 10.^[2.7 1.75 2.9 2.5];
%modelt = [20000 30000 350000];
%mmodel = [modelr modelt];

modelr = [500 500 500 500]; % initial resistivity (Ohm-m)
modelt = [20000 20000 20000]; % iniial thickness (m)
mmodel = [modelr modelt];

%logFrequencies = -3:0.2:3;
%frequency = 10.^logFrequencies;

m_est = [2.5970e+02   8.4931e+01   2.9850e+02   3.4647e+02   2.0151e+04   2.0026e+04   1.9877e+04];
m = m_est;

lr = length(modelr);
lt = length(modelt);

r_est = m_est(1:lr);
t_est = m_est(1+lr:lr+lt);

[rhoa_calc, phase_calc] = modelagem1DMT(r_est, t_est, frequency)

subplot(2,2,1)

loglog(frequency,rhoa_calc,'-','color','b','LineWidth',2)
hold on
loglog(frequency,rhoa, '.','color','r','MarkerSize',15)
grid on
axis tight
xlabel('\bf \fontsize{10}\fontname{Times}Frequência (Hz)');
ylabel('\bf \fontsize{10}\fontname{Times}Resistividade Aparente (Ohm.m)');
%title(['\bf \fontsize{12}\fontname{Times}RESPONS - LM || RMS : ', num2str(rms_error),' || iter : ', num2str(iter)]);
leg = legend('Resistividade Aparente calculada','Resistividade Aparente observada');
%set(leg,'Location','South','fontsize',8);
set(leg,'fontsize',8);
hold off
subplot(2,2,2)
loglog(frequency,phase_calc,'-','color','b','LineWidth',2)
hold on
loglog(frequency,faserealradianos, '.','color','r','MarkerSize',15)
grid on
xlabel('\bf \fontsize{10}\fontname{Times}Frequência (Hz)');
ylabel('\bf \fontsize{10}\fontname{Times}Fase (rad)');
%title(['\bf \fontsize{12}\fontname{Times}RESPONS - LM || RMS : ', num2str(rms_error),' || iter : ', num2str(iter)]);
leg = legend('Fase calculada ','Fase observada');
axis([0.0004 0.0400]);
%set(leg,'Location','South','fontsize',8);
set(leg,'fontsize',8);
hold off

r = m(1:lr);
t = m(1+lr:lr+lt);
rr = [0,r];
tt = [0,cumsum(t),max(t)*5];
modelrr = [0,modelr];
modeltt = [0,cumsum(modelt),max(modelt)*5];

subplot(2,2,3);
hold off
stairs(modelrr,modeltt,'-','color','r','LineWidth',2);
hold on
stairs(rr,tt,'b','LineWidth',2);
% '.','color','r','MarkerSize',15
set(gca,'Ydir','reverse');
set(gca,'Xscale','log');
leg = legend('Inicial','Estimado');
%set(leg,'Location','South','fontsize',8);
set(leg,'fontsize',8);
%axis([1 100000 0 100000]);
grid on
xlabel('Resistividade(Ohm-m)','fontweight','bold','fontsize',10);
ylabel('Profundidade (m)','fontweight','bold','fontsize',10);
%set(gca,'XTick',[1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4]);

subplot(2,2,4);
plot(1:15,rms_error(1:15).^2,'-','color','b','LineWidth',2)
xlabel('Iteração','fontweight','bold','fontsize',10);
ylabel('Função objetivo','fontweight','bold','fontsize',10);
axis([1 16]);
grid on






