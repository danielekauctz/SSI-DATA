%% SSI-DATA APPLICATION

% Numerical implementation: 
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

% References:
% CARINI, M. R. (2021). Identificação das Propriedades Dinâmicas de Estruturas submetidas a Ações Ambientais. 
  ...Doctoral Thesis, Federal University of Rio Grande do Sul, Brazil.
% PEETERS, B. (2000). System Identification and Damage Detection in Civil Engineering. 
  ...PhD thesis, Katholieke Universiteit Leuven, Belgium.
      
% Parameters:
% y: output matrix 
% i: half the number of block rows of the output data Hankel matrix (i x l > no)
% no: state-space model order
% sf: sampling frequency
% FN, Phi, ZT: natural frequencies, mode shapes and damping ratios of the system

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
tol = [0.01 0.01 0.10];    % [fn zt phi] 
x_lim = [0.01 10];         % frequency range

[stbf,stbz,stbm] = stabilization(FN,ZT,PHI,tol,ordem,x_lim);  
