function [EF] = achtpunktalgorithmus(Korrespondenzen, K)
    % Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
    % mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
    % vorliegt oder nicht
    
    %% (3.1)
    % Transformieren sie die Korrespondenzpunktpaare in 2 Arrays mit den homogenen kalibirierten Koordinaten x1 und x2
    N = size(Korrespondenzen, 2);
    
    x1 = [Korrespondenzen(1:2,:);
          ones(1,N)];
    x2 = [Korrespondenzen(3:4,:);
          ones(1,N)];
    if exist('K','var')
        % Kalibrierungsmatrix K liegt vor, Koordinaten umrechnen
        %k1 = inv(K);
        x1 = K\x1;
        x2 = K\x2;
    end
    
    % A ist the Kronecker Product of the coordinates vectors x1 and x2
    A = zeros(N,9);
    for i=1:N
        A(i,:) = kron(x1(:,i), x2(:,i))';
    end
    %size(A)
    
    % Berechne SVD von A
    [~,~,V] = svd(A);
    
    %% Schaetzung der Matrizen (3.2)
    % 9. Spalte von V der SVD minimiert das Problem: argmin_(||E_s||2=1) || AE_s ||_2^2
    Gs = V(:,9);
    
    % Umsortieren zu 3x3 Matrix
    G = reshape(Gs,3,3);
    
    % Finde die n�chste Essentielle Matrix zu G: argmin_(E) || E-G ||_2^2
    % Singulaerwertzerlegung der umsortierten Loesung
    [Ug, Sg, Vg] = svd(G);
    
    % Zu-Eins-Setzen der ersten beiden Singulaerwerte und zu null-Setzen des kleines Singulaerwerts
    %sigma = (Sg(1,1) * Sg(2,2)) / 2;    
        
    if exist('K','var')
        % Kalibrierungsmatrix K liegt vor, Essentielle Matrix zurueck geben
        S = diag([1,1,0]);
        F = Ug * S * Vg';
        EF = F;
    else
        % Kalibirierungsmatrix K liegt nicht vor, Fundamentalmatrix zurueck geben
        S = Sg;
        S(3,3) = 0;
        F = Ug * S * Vg';
        EF = F;
    end
end