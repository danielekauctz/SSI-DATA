%% Data-driven stochastic subspace identification (SSI-DATA)

% Numerical implementation: 
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

% References:
% FADEL MIGUEL, L. F. (2007). Identificaçao de sistemas e avaliaçao da integridade de estruturas treliçadas. 
  ...Doctoral Thesis, Federal University of Rio Grande do Sul, Brazil.
% PEETERS, B. (2000). System Identification and Damage Detection in Civil Engineering. 
  ...PhD thesis, Katholieke Universiteit Leuven, Belgium.
% VAN OVERSCHEE, P.; DE MOOR, B. (1996). Subspace identification for linear systems: Theory - Implementation - Applications. 
  ...Leuven: Kluwer Academic Publishers, 254 p.

% Parameters:
% y: output matrix 
% i: half the number of block rows of the output data Hankel matrix (i x l > no)
% no: state-space model order
% sf: sampling frequency
% dt: sampling period
% l: number of outputs
% n: number of time samples
% N: number of columns in the Hankel matrix
% H: Hankel matrix
% Oi: Observability matrix
% Pi: Projection matrix
% C: Output matrix
% A: State matrix
% FN, Phi, ZT: natural frequencies, mode shapes and damping ratios of the system

function [FN,PHI, ZT] = SSI_DATA(y,i,no,sf)

dt = 1/sf;

[l,n] = size(y);
j = n-2*i+1;
yy = y/sqrt(j);

clear y

% Hankel matrix
H = zeros(l*2*i,j);
for k = 1:2*i
	H((k-1)*l+1:k*l,:) = yy(:,k:k+j-1);
end

clear yy

% QR factorization and projection
R = triu(qr(H'))';
R = R(1:2*i*l,1:2*i*l);

clear H

Pi = R(l*i+1:2*l*i,1:l*i);

clear R

% Singular value decomposition (SVD) 
[U,S,~] = svd(Pi);
ss = diag(S);

clear Pi 
clear S

% Modal properties for different model order
for ii = 1:no

  U1 = U(:,1:ii);
  Oi = U1*diag(sqrt(ss(1:ii)));
  Oi1 = U1(1:l*(i-1),:)*diag(sqrt(ss(1:ii)));
  
  % Matrices A & C
  A = pinv(Oi1)*Oi(l+1:l*i,:); 
  C = Oi(1:l,:);
  
  clear U1
  clear Oi
  clear Oi1
  
  [Psi,Lambda] = eig(A);  
  fn  = zeros(1,length(Lambda));
  csi = zeros(1,length(Lambda));
  phi  = C*Psi;
  
   for jj = 1:length(Lambda)
       u(jj) = log(Lambda(jj,jj))/dt;
       fn(jj) = abs(u(jj))/(2*pi);
       csi(jj) = -real(u(jj))/abs(u(jj));
   end
   
  FD  = [fn(:) csi(:)];
  [~,idx] = sort(FD(:,1));
  FDS = FD(idx,:);
  ww=find(FDS(:,1)>0 & FDS(:,2)>0);
  FDTf(:,ii)=[FDS(ww,1);zeros(no-length(ww),1)];
  FDTz(:,ii)=[FDS(ww,2);zeros(no-length(ww),1)];
  
  [~,cc] = size(phi);
   for kk = 1:cc
       ms(:,kk) = (abs(phi(:,kk))./max(abs(phi(:,kk)))).*sign(real(phi(:,kk)));
   end
   ms1 = ms(:,idx);
   ms2 = ms1(:,ww);
   [ee,ff] = size(ms2);
   FM1((ee*ii)-(ee-1):ee*ii,:) = [ms2 zeros(ee,no-ff)]; 

end

FN  = FDTf;
PHI  = FM1;
ZT = FDTz;

end