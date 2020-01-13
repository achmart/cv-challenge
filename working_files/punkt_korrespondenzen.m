function Korrespondenzen = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,winkel_1,winkel_2,do_rotate,varargin)
    % In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
    % Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
    % zu ermitteln.
    
    %% Input parser (2.1)
    p = inputParser;
    addRequired(p,'I1');
    addRequired(p,'I2');
    
    default_window = 25;
    default_min_corr = 0.95;
    default_plot = false;
    
    check_window = @(x) isscalar(x) && isnumeric(x) && (x > 1) && (mod(x,2) == 1);
    check_min_corr = @(x) isscalar(x) && isnumeric(x) && (x > 0.0) && (x < 1.0);
    check_do_plot = @(x) isscalar(x) && islogical(x);
    
    addParameter(p, 'window_length', default_window, check_window); 
    addParameter(p, 'min_corr', default_min_corr, check_min_corr);
    addParameter(p, 'do_plot', default_plot, check_do_plot);
    
    parse(p, I1, I2, varargin{:});
    
    %Im1 = double(p.Results.I1);
    %Im2 = double(p.Results.I2);
    
    %% Save variables
    window_length = p.Results.window_length;
    min_corr = p.Results.min_corr;
    do_plot = p.Results.do_plot;
    
    %% Merkmalsvorbereitung (2.2)
    [H1,W1] = size(I1);
    w = floor(ceil(sqrt(2)*window_length) / 2) + 1;
    
    A = floor(w);
    C = floor(w);
    B = ceil(W1 - w);
    D = ceil(H1 - w);
    
    no_pts1 = size(Mpt1,2);
    no_pts2 = size(Mpt2,2);
    
    % Gehe durch alle Merkmale von Bild I1
    maxIndex = no_pts1;
    i = 1;
    while(i <= maxIndex)
        % Hole Koordinaten des Merkmals
        XY = Mpt1(:,i);
        x = XY(1);
        y = XY(2);
        % �berpr�fe die Distanz zum Rand
        if(x <= A || x >= B || y <= C || y >= D)
            % Entferne dieses Merkmal
            Mpt1(:,i) = [];
            no_pts1 = no_pts1 - 1;
            % Don't increment
        else
            %fprintf('Not deleting point x: %d and y : %d\n', x, y);
            % Increment only here
            i = i + 1;
        end
        maxIndex = size(Mpt1,2);
    end
    
    % Gehe durch alle Merkmale von Bild I2
    maxIndex = no_pts2;
    i = 1;
    while(i <= maxIndex)
        % Hole Koordinaten des Merkmals
        XY = Mpt2(:,i);
        x = XY(1);
        y = XY(2);
        % �berpr�fe die Distanz zum Rand
        if(x <= A || x >= B || y <= C || y >= D)
            % Entferne dieses Merkmal
            Mpt2(:,i) = [];
            no_pts2 = no_pts2 - 1;
            % Don't increment
        else
            %fprintf('Not deleting point x: %d and y : %d\n', x, y);
            % Increment only here
            i = i + 1;
        end
        maxIndex = size(Mpt2,2);
    end
    
    
    %% Normierung der Fenster (2.3)
    Mat_feat_1 = [];
    Mat_feat_2 = [];
    %Mat_feat_1 = zeros(size(Mpt1,2));
    
    r = floor(window_length / 2);
    % New patch edge length
    % (to include the whole window in the rotated patch)
    diag = ceil(sqrt(2)*window_length);
    r2 = floor(diag / 2);
    cx = r2;
    cy = r2;
    
    for i=1:size(Mpt1,2)
        XY = Mpt1(:,i);
        x = XY(2);
        y = XY(1);
        if(do_rotate)
            A = I1((x-r2):(x+r2), (y-r2):(y+r2));
            A = double(A);
            % Rotieren um den Winkel
            A = imgrotate(A,-winkel_1(x,y));
            %size(A)
            [cy2, cx2] = size(A);
            cy = ceil(cy2 / 2);
            cx = ceil(cx2 / 2);
            A = A( (cy-r):(cy+r), (cx-r):(cx+r) );
        else
            A = I1((x-r):(x+r), (y-r):(y+r));
            A = double(A);
        end
        A_mean = mean(mean(A));
        A_std = sqrt( (1/(window_length^2-1)) * trace((A-A_mean)'*(A-A_mean)) );
        Anorm = (A - A_mean) / A_std;
        if(isnan(Anorm))
            display([x, y]);
            Anorm = zeros(size(Anorm));
        end
        Mat_feat_1 = [Mat_feat_1, Anorm(:)];
    end
    
    for i=1:size(Mpt2,2)
        XY = Mpt2(:,i);
        x = XY(2);
        y = XY(1);
        if(do_rotate)
            A = I2((x-r2):(x+r2), (y-r2):(y+r2));
            A = double(A);
            % Rotieren um den Winkel
            A = imgrotate(A,-winkel_1(x,y));
            %size(A)
            [cy2, cx2] = size(A);
            cy = ceil(cy2 / 2);
            cx = ceil(cx2 / 2);
            A = A( (cy-r):(cy+r), (cx-r):(cx+r) );
        else
            A = I2((x-r):(x+r), (y-r):(y+r));
            A = double(A);
        end
        A_mean = mean(mean(A));
        A_std = sqrt( (1/(window_length^2-1)) * trace((A-A_mean)'*(A-A_mean)) );
        Anorm = (A - A_mean) / A_std;
        if(isnan(Anorm))
            display([x, y]);
            Anorm = zeros(size(Anorm));
        end
        Mat_feat_2 = [Mat_feat_2, Anorm(:)];
    end
    
    %% NCC Brechnung (2.4)
    S1 = size(Mat_feat_1, 2);
    S2 = size(Mat_feat_2, 2);
    NCC_matrix = zeros(S2, S1);
    N = window_length^2;
    f = 1/(N-1);
    
    for i=1:S2
        for j=1:S1
            %INDEX = [i,j]
            WnVn = f * Mat_feat_1(:,j)' * Mat_feat_2(:,i);
            if(isnan(WnVn))
                disp('error" WnVn is NaN, (i,j)=...')
                disp([i,j])
                NCC_matrix(i,j) = 0.0;
            end
            NCC_matrix(i,j) = WnVn;    
        end
    end
    
    %NCC_matrix(200:end,:)
    
    % Alle Eintraege zu 0 setzen, die unter min_corr liegen
    NCC_matrix(NCC_matrix < min_corr) = 0;
    
    % Sortieren aller Korrespondenzwerte, die ungleich 0 sind in absteigender Reihenfolge der Korrespondenzstaerke
    [RS,CS,v] = find(NCC_matrix);
    linearInd = sub2ind(size(NCC_matrix), RS, CS);
    [~, I] = sort(v, 'descend');
    sorted_index = linearInd(I);
    
    %% Korrespondenzen (2.5)
    Korrespondenzen = [];
    %[S2, S1] = size(NCC_matrix) % S2 (S1): features in Mat_feat_2 (1)
    size(NCC_matrix);
    i = 1;
    maxLength = size(sorted_index,1);
    [H,L] = size(I1);
    
    while(i <= maxLength)
        idx = sorted_index(i);
        % Ueberpruefe ob dieser Wert null ist
        [x,y] = ind2sub(size(NCC_matrix), idx);
        x_dist = abs(Mpt1(1,y)-Mpt2(1,x));
        y_dist = abs(Mpt1(2,y)-Mpt2(2,x));
        
        if((NCC_matrix(x,y) ~= 0) && (x_dist < L/6) &&(y_dist<0.1*H))
            %col1 = Mat_feat_1(:,y);
            %col2 = Mat_feat_2(:,x);
            
            % Speichere diese x,y Koordinaten
            vec = [Mpt1(1,y);
                   Mpt1(2,y);
                   Mpt2(1,x);
                   Mpt2(2,x)];
            
            % Anhaengen an Gesamt-Matrix
            Korrespondenzen = [Korrespondenzen, vec];
            
            % Spalte in NCC_matrix gleich 0 setzen
            NCC_matrix(:,y) = 0;
            NCC_matrix(x,:) = 0;
        end
        i = i + 1;
    end
    
    
    
    %% Zeige die Korrespondenzpunktpaare an (2.6)
    if(do_plot)
        figure
        imshow(I1)
        alpha 0.5
        hold on
        imshow(I2)
        alpha 0.5
        plot(Korrespondenzen(1,:), Korrespondenzen(2,:),'or')
        plot(Korrespondenzen(3,:), Korrespondenzen(4,:),'xg')
        for i=1:size(Korrespondenzen,2)
            plot(Korrespondenzen([1, 3],i), Korrespondenzen([2,4],i),'-')
        end
        hold off
        title('All correspondences');
    end
end
