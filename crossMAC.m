%% MODAL ASSURANCE CRITERION (MAC)

% Numerical implementation: 
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

% Reference:
% ALLEMANG, R. J.; BROWN, D. L. (1982). A Correlation Coefficient for Modal Vector Analysis. 
  ...Proceedings of the 1st International Modal Analysis Conference (IMAC), USA, p. 110-116.

% Parameters:
% Phi1, Phi2: mode shapes matrices
        ... rows = degrees of freedom
        ... columns = number of modes
% n1, n2: number of modes of each matrix
% plot: if plot = 1 generates bar chart of MAC

function MAC = crossMAC(Phi1,Phi2,plot)

    [~, n1]= size(Phi1);
    [~, n2]= size(Phi2);

    MAC = zeros(n1,n2,'double');

    for mode1 = 1:n1    
        for mode2 = 1:n2
            MAC(mode1,mode2) = abs(Phi1(:,mode1)'*Phi2(:,mode2))/ ...
                               sqrt(abs((Phi1(:,mode1)'*Phi1(:,mode1))*(Phi2(:,mode2)'*Phi2(:,mode2))));                          
        end   
    end
    
    if plot == 1
        figure
        barMAC = bar3(MAC);

        colormap(jet);
        title('Modal Assurance Criterion','fontsize',15)
    end

end