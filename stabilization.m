%% STABILIZATION OF NATURAL FREQUENCIES, MODE SHAPES AND DAMPING RATIOS

% Numerical implementation: 
% Gustavo Zeni & Daniele Kauctz Monteiro (2022)

% References:
% CARINI, M. R. (2021). Identificação das Propriedades Dinâmicas de Estruturas submetidas a Ações Ambientais. 
  ...Doctoral Thesis, Federal University of Rio Grande do Sul, Brazil.
% PEETERS, B. (2000). System Identification and Damage Detection in Civil Engineering. 
  ...PhD thesis, Katholieke Universiteit Leuven, Belgium.

% Parameters:
% stbf: position of stable frequencies
% stbz: position of stable damping ratios
% stbm: position of stable mode shapes
  
function [stbf,stbz,stbm] = stabilization(FN,ZT,PHI,tol,ordem,x_lim)

P(:,:,1) = FN;
P(:,:,2) = ZT;

for i = 1:2
    out(:,:,i) = stab(P(:,:,i),tol(i));
end

stbf = squeeze(out(:,:,1));
stbz = squeeze(out(:,:,2));
stbm = stab_modes(PHI,tol(3));

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
xlim(x_lim)

end

function out = stab(F,tol)

F = squeeze(F(:,:,1));
stb = zeros(size(F,1),size(F,2));
for ord = 1:size(F,2)-1
    [row,~,v] = find(F(:,ord+1));
    if ~isempty(v)
        for ii = 1:length(v)
            lim = v(ii)*[(1-tol) (1+tol)];
            [row2,~] = find(F(:,ord)>=lim(1) & F(:,ord)<= lim(2));
            if ~isempty(row2)
                stb(row(ii),ord+1) = 1;
            else
                %do nothing
            end
        end
    else
        %do nothing
    end
end      
out = stb;
end

function out = stab_modes(F,tol)

n = size(F,1)/size(F,2);
stb = zeros(size(F,2));

for ord = 1:size(F,2)-1
    phi1 = F(n*ord+1:(ord+1)*n,:);
    phi2 = F(n*(ord-1)+1:ord*n,:);
    if ~isempty(phi1) && ~isempty(phi2)
        mac = MAC(phi1,phi2);
        mac(isnan(mac)) = 0;
        [~,col] = find(mac>=1-tol);
        stb(col,ord+1) = 1;             
    end
                
end

out = stb;

end
