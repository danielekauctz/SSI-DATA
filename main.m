%% SSI-DATA APPLICATION

% Numerical implementation: 
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

% References:
% CARINI, M. R. (2021). Identificação das Propriedades Dinâmicas de Estruturas submetidas a Ações Ambientais. 
  ...Doctoral Thesis, Federal University of Rio Grande do Sul, Brazil.
% PEETERS, B. (2000). System Identification and Damage Detection in Civil Engineering. 
  ...PhD thesis, Katholieke Universiteit Leuven, Belgium.

%% Carini's Acceleration signal
yk = importdata('yk.txt'); 

%% SSI-DATA
% i x l > no
no = 20;
i = 10;
sf = 25;

[FN,PHI,ZT] = SSI_DATA(yk,i,no,sf);
ordem = 1:size(FN,2);

%% Stabilization Diagram
close all
tol = [0.01 0.01 0.10]; % [fn zt phi] 
[stbf,stbz,stbm] = stabilization(FN,ZT,PHI,tol);

F1 = (FN.*stbf)';              % stable frequencies (red)
F2 = (FN.*stbf.*stbm)';        % stable frequencies & mode shapes (blue)
F3 = (FN.*stbf.*stbm.*stbz)';  % stable frequencies, mode shapes & damping ratios (green)

figure
plot(FN',ordem,'k.');hold on
plot(F1,ordem,'ro')
plot(F2,ordem,'bo')
plot(F3,ordem,'go');hold off
xlabel('Frequency [Hz]','FontSize',12);ylabel('Model order','FontSize',12)
title('Stabilization Diagram','FontSize',14)
grid on
xlim([0.01 10])

figure
plot(ZT',ordem,'k.');hold on
plot((ZT.*stbz)',ordem,'ro');hold off
xlabel('Damping ratio','FontSize',12);ylabel('Model order','FontSize',12)
title('Stabilization Diagram','FontSize',14)
grid on
xlim([0.0001 0.1])
