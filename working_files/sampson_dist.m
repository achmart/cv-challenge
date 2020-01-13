function sd = sampson_dist(F, x1_pixel, x2_pixel)
    % Diese Funktion berechnet die Sampson Distanz basierend auf der
    % Fundamentalmatrix F
    
    %% Sampson-Distanz (3.5)
    e3_dach = dach([0,0,1]);
    Fx1 = F * x1_pixel;
    zaehler = diag(x2_pixel' * Fx1);
    zaehler = (zaehler.^2)';
    %sizeZaehler = size(zaehler);    
    nenner1 = vecnorm(e3_dach * F * x1_pixel);
    nenner2 = vecnorm((x2_pixel' * F * e3_dach)');
    nenner = nenner1.^2 + nenner2.^2;
    %sizeNenner = size(nenner);    
    sd = zaehler ./ nenner;
end