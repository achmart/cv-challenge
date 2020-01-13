function [T1, R1, T2, R2, U, V] = TR_aus_E(E)
    % Diese Funktion berechnet die moeglichen Werte fuer T und R
    % aus der Essentiellen Matrix
    
    % Berechnung der SVD von E (4.1)
    [U,S,V] = svd(E);
    
    % Stellen sie sicher, dass U und V Rotationsmatrizen sind
    % Hinreichende Eigenschaft: Quadritsche Matrix mit zu 1 normierte, orthogonoalen Spalten, Det(A)=1
    detU = det(U);
    detV = det(V);
    
    if(detU < 0)
        U = U * diag([1, 1, -1]);
    end
    
    if(detV < 0)
        V = V * diag([1, 1, -1]);
    end
    
    % Normierung der Spalten
    nU = vecnorm(U);
    nV = vecnorm(V);
    
    U = U ./ nU;
    V = V ./ nV;
    
    %detU = det(U);
    %detV = det(V);
    
    Rz1 = [0 -1 0; +1 0 0; 0 0 1];
    Rz2 = [0 +1 0; -1 0 0; 0 0 1];
    
    R1 = U * Rz1' * V';
    R2 = U * Rz2' * V';
    
    T1_dach = U * Rz1 * S * U';
    T2_dach = U * Rz2 * S * U';
    T1 = [T1_dach(3,2); T1_dach(1,3); T1_dach(2,1)];
    T2 = [T2_dach(3,2); T2_dach(1,3); T2_dach(2,1)];
    %[T1,R1,T2,R2,U,V]; 
end