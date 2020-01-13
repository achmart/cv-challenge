function [T, R, lambda, P1, camC1, camC2] = rekonstruktion(T1, T2, R1, R2, Korrespondenzen, K, do_plot)

    %% Preparation (4.2)
    N = size(Korrespondenzen, 2);
    %K_inv = inv(K);
    x1 = K\[Korrespondenzen(1:2,:); ones(1,N)];
    x2 = K\[Korrespondenzen(3:4,:); ones(1,N)];
    
    T_cell = { T1, T2, T1, T2 };
    
    R_cell = { R1, R1, R2, R2 };
    
    % Cell array d_cell definiert mit { ... }
    d_cell = {zeros(N,2), zeros(N,2), zeros(N,2), zeros(N,2)};
    
    %% Rekonstruktion Gleichungssysteme (4.3)
    
    N = size(Korrespondenzen, 2);
    
    % Iterieren Sie druch alle Kombinationen
    for i=1:4 
        %ITERATION = i
        % Matrizen M1 und M2 erstellen
        T = cell2mat(T_cell(i));
        R = cell2mat(R_cell(i));
        
        M1 = zeros(3*N, N+1);
        M2 = zeros(3*N, N+1);
        %sizeM1 = size(M1);
        %sizeM2 = size(M2);
        
        for k=1:N
           % M1
           M1( (3*k-2):(3*k), k) = dach(x2(:,k)) * R * x1(:,k);
           M1( (3*k-2):(3*k), N+1) = dach(x2(:,k)) * T;
           % M2           
           M2( (3*k-2):(3*k), k) = dach(x1(:,k)) * R' * x2(:,k);
           M2( (3*k-2):(3*k), N+1) = -dach(x1(:,k)) *R' *T;
        end       
        
        % Bestimmen sie die Loesung der Gleichungssysteme mittels SVD
        [U1,S1,V1] = svd(M1);
        [U2,S2,V2] = svd(M2);
        
        % Letzte Spalte von V1 / V2 loest das Minimierungsproblem min ||M*d|| statt M*d=0
        % Loesung des Minimierungsproblems ist die letzte Spalte von V1 / V2
        d1 = V1(:,end);
        d2 = V2(:,end);
        
        % Normierung sodass der letzte Wert 1 ergibt
        d1 = d1 ./ d1(end);
        d2 = d2 ./ d2(end);
        
        % Kopieren der Werte (bis auf den letzten Wert) in d_cell array
        d = [d1(1:end-1), d2(1:end-1)];
        %SIZE_D = size(d)        
        d_cell(i) = mat2cell(d, size(d,1), size(d,2));
        % Nur zum test ob das speichern in cell-Array geklappt hat
        %d_test = cell2mat(d_cell(i))
    end
    
    % Finde Index i, fuer den die meisten positiven Tiefeninformationen in d_cell(i) enthalten sind
    d_nr_pos = [];
    for j=1:4
       d_temp = cell2mat(d_cell(j));
       nr_pos = (d_temp >= 0);
       % reduce to 1 dimension
       nr_pos = nr_pos(:);       
       d_nr_pos = [d_nr_pos, sum(nr_pos)];
    end
           
    [~, i_max] = max(d_nr_pos);
    
    T = cell2mat(T_cell(i_max));
    R = cell2mat(R_cell(i_max));
    lambda = cell2mat(d_cell(i_max));
    
    %sizeLambda = size(lambda);
    %sizeM1 = size(M1);
    %sizeM2 = size(M2);
    %{T, R, lambda, M1, M2}
    
    %% (4.4)
    
    % Berechnen Sie die Weltkoordinaten der Bildpunkte P1 aus den Tielfeninformationen lambda und den Bildkoordinaten x1.
    size(lambda);
    size(x1);
    N = size(x1,2);

    P1 = lambda(:,1)' .* x1;
    % Ecken des Kameraframe1
    camC1 = [-0.2 0.2 0.2 -0.2; 0.2, 0.2, -0.2, -0.2; 1 1 1 1];

    % Berechnen sie mithilfe dieser Matrix und der Euklidischen Bewegung die Weltkoordinaten des Cameraframes 2
    camC2 = R' * (camC1 - T);
    
        % Zeichnen Sie die Weltkoordinaten P1 mit Hilfe der Funktion scatter3()
    if do_plot == true
        figure
        scatter3(P1(1,:), P1(2,:), P1(3,:))
        grid on
        hold on
        axis equal
        xlabel('X');
        ylabel('Y');
        zlabel('Z');

        % und numerieren Sie diese im Bild mit der text() Funktion.
        indices = cellstr(num2str([1:N]'));
        text(P1(1,:), P1(2,:), P1(3,:), indices);



        % Zeichnen von camC1 in blau
        plot3(camC1(1,:), camC1(2,:), camC1(3,:), 'b');
        text(camC1(1,1), camC1(2,1), camC1(3,1), 'CamC1');

        % Zeichnen von camC2 in rot
        plot3(camC2(1,:), camC2(2,:), camC2(3,:), 'r');
        text(camC2(1,1), camC2(2,1), camC2(3,1), 'CamC2');
        hold off

        % Kamerposition und Neigung
        campos([43,-22,-87]);
        camup([0, -1, 0]);
    end
    
    %{T, R, lambda, P1, camC1, camC2};
end

