function [repro_error, x2_repro] = rueckprojektion(Korrespondenzen, P1, Image2, T, R, K)
    % Diese Funktion berechnet den mittleren Rueckprojektionsfehler der 
    % Weltkooridnaten P1 aus Bild 1 im Cameraframe 2 und stellt die 
    % korrekten Merkmalskoordinaten sowie die rueckprojezierten grafisch dar.
    % Bestimmen sie zuerst die Weltkoordinaten von P1 im Cameraframe2
    
    %% (4.5)
    
    P2 = R * P1 + T;
    
    %  Berechnen sie die homogenen Bildkoordinaten x2_repro
    x2_repro = K * P2;
    % Normieren auf Z=1
    x2_repro = x2_repro ./ x2_repro(3,:);
    
    x2 = Korrespondenzen(3:4,:);
    N = size(x2,2);
    
    % Zeigen sie bild 2 an und zeichnen sie die gefundenen Merkmale sowie x2repro ein
    figure
    imshow(Image2)
    hold on
    scatter(x2_repro(1,:), x2_repro(2,:), 'or')
    indices = cellstr(num2str((1:N)'));
    text(x2_repro(1,:), x2_repro(2,:), indices, 'Color','red')
    
    scatter(x2(1,:), x2(2,:), 'xg')
    text(x2(1,:), x2(2,:), indices, 'Color','green')
    for i=1:size(x2,2)
        plot( [x2_repro(1,i), x2(1,i)], [x2_repro(2,i), x2(2,i)],'-y')
    end
    hold off
        
    repro_error = 1/N * sum(vecnorm([x2; ones(1,N)] - x2_repro));
    
    %{repro_error, x2_repro};
end

