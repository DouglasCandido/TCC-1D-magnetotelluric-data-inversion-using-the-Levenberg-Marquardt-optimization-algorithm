close all;

clc;

clear;

load rms_error.mat;
load rhoa_ruido2.mat;
load phase_ruido2.mat;

% Load model for forward modeling
%modelr = [30000 5000 1000]; % modelo de resistividade verdadeiro (Ohm-m)
%modelt = [15 18]; % modelo de espessura verdadeiro (Km)
modelr = [1000 100 1000]; % modelo de resistividade verdadeiro (Ohm-m)
modelt = [1000 1000]; % modelo de espessura verdadeiro (m)
mmodel = [modelr modelt]; % modelo verdadeiro

logFrequencies = -3:0.2:3;
frequency = 10.^logFrequencies;

m_est = [1000.840     74.298    971.515   1066.646    686.065];
m = m_est;

lr = length(modelr);
lt = length(modelt);

r_est = m_est(1:lr);
t_est = m_est(1+lr:lr+lt);

[rhoa_calc, phase_calc] = modelagem1DMT(r_est, t_est, frequency);

subplot(2,2,1)

loglog(frequency,rhoa_calc,'-','color','b','LineWidth',2)
hold on
loglog(frequency,rhoa_ruido2, '.','color','r','MarkerSize',15)
grid on
axis tight
xlabel('\bf \fontsize{10}\fontname{Times}Frequência (Hz)');
ylabel('\bf \fontsize{10}\fontname{Times}Resistividade Aparente (Ohm.m)');
%title(['\bf \fontsize{12}\fontname{Times}RESPONS - LM || RMS : ', num2str(rms_error),' || iter : ', num2str(iter)]);
leg = legend('Resistividade Aparente calculada','Resistividade Aparente Observada');
%set(leg,'Location','South','fontsize',8);
set(leg,'fontsize',8);
hold off
subplot(2,2,2)
loglog(frequency,phase_calc,'-','color','b','LineWidth',2)
hold on
loglog(frequency,phase_ruido2, '.','color','r','MarkerSize',15)
grid on
xlabel('\bf \fontsize{10}\fontname{Times}Frequência (Hz)');
ylabel('\bf \fontsize{10}\fontname{Times}Fase (rad)');
%title(['\bf \fontsize{12}\fontname{Times}RESPONS - LM || RMS : ', num2str(rms_error),' || iter : ', num2str(iter)]);
leg = legend('Fase calculada ','Fase Observada');
%set(leg,'Location','South','fontsize',8);
set(leg,'fontsize',8);

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
leg = legend('Verdadeiro','Estimado');
%set(leg,'Location','South','fontsize',8);
set(leg,'fontsize',8);
%axis([1 15000 0 2000]);
grid on
xlabel('Resistividade (Ohm-m)','fontweight','bold','fontsize',10);
ylabel('Profundidade (m)','fontweight','bold','fontsize',10);
set(gca,'XTick',[1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4]);

subplot(2,2,4);
plot(1:10,rms_error.^2,'-','color','b','LineWidth',2)
xlabel('Iteração','fontweight','bold','fontsize',10);
ylabel('Função objetivo','fontweight','bold','fontsize',10);
grid on








