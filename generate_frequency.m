clc;

clear;

x = [0.001122018	2137.96209	10
0.002924152	1200	10
0.007943282	1000	10
0.015625	850	10
0.035645113	700	10
0.061376201	389.045145	10
0.159955803	177.827941	10
0.534564359	64.5654229	10
4.073802778	19.498446	10
14.99684836	42	10
44.44444444	105	10
100	323.5936569	10
301.995172	489.7788194	10
1472.312502	316.227766	10
4444.444444	250	10];

frequency = x(:,1);

save("frequency.mat", "frequency")