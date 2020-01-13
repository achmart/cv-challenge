function [merkmale,winkel] = harris_detektor(input_image, varargin)
    % In dieser Funktion soll der Harris-Detektor implementiert werden, der
    % Merkmalspunkte aus dem Bild extrahiert

    %% Input parser (Upgraded 1.7)
    p = inputParser;
    addRequired(p,'input_image');
    
    default_segment = 15;
    default_k = 0.05;
    default_tau = 1e6;
    default_plot = false;
    default_min_dist = 20;
    default_tile_size = [200, 200];
    default_N = 5;

    check_segment = @(x) isscalar(x) && isnumeric(x) && (x > 1) && (mod(x,2) == 1);
    check_k = @(x) isscalar(x) && isnumeric(x) && (x >= 0.0) && (x <= 1.0);
    check_tau = @(x) isscalar(x) && isnumeric(x) && (x > 0);
    check_do_plot = @(x) isscalar(x) && islogical(x);
    check_min_dist = @(x) isscalar(x) && isnumeric(x) && (x >= 1);
    check_tile_size = @(x) isnumeric(x);
    check_N = @(x) isscalar(x) && (x >= 1);

    addParameter(p, 'segment_length', default_segment, check_segment)
    addParameter(p, 'k', default_k, check_k)
    addParameter(p, 'tau', default_tau, check_tau)
    addParameter(p, 'do_plot', default_plot, check_do_plot)
    addParameter(p, 'min_dist', default_min_dist, check_min_dist)
    addParameter(p, 'tile_size', default_tile_size, check_tile_size)
    addParameter(p, 'N', default_N, check_N)
    
    parse(p, input_image, varargin{:})
    
    if(isscalar(p.Results.tile_size))
        tile_size = [p.Results.tile_size, p.Results.tile_size];
    else
        tile_size = p.Results.tile_size;
    end
    
    %% Assign Parser results to corresponding variables
    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    min_dist = p.Results.min_dist;
    %N = p.Results.N;  
    
    
    %% Vorbereitung zur Feature Detektion (1.4)
    % Pruefe ob es sich um ein Grauwertbild handelt
    %dims = size(input_image);
    if ndims(input_image) > 2
        error("Image format has to be NxMx1");
    end
    
    % Konvertiere zu double
    % Rette img
    %img = input_image;
    %input_image = double(input_image);

    % Approximation des Bildgradienten
    [Ix, Iy] = sobel_xy(double(input_image));
    
    % Gewichtung
    % Gaussfilter
    N = segment_length;
    indices = (1:segment_length)';
    indices = indices - ceil(N/2);
    % Gewichtungs-Vektor w
    w = exp(-0.5*(2.5*indices / (N/2) ).^2);
    % Normierungskonstante C
    C = 1 / sum(w);
    w = C * w;
    W = w*w';

    % Harris Matrix G    
    G11 = conv2(Ix.^2, W, 'same');
    G12 = conv2(Ix.*Iy, W, 'same');
    G22 = conv2(Iy.^2, W, 'same');
    
    % Ix, Iy filtern
    Ix_filt = conv2(Ix, W, 'same');
    Iy_filt = conv2(Iy, W, 'same');
    
    %% Merkmalsextraktion ueber die Harrismessung
    [M, N] = size(input_image);
    H = zeros(M, N);
    winkel = zeros(M, N);
    
    % This code has been vectorized below
    %for i=1:M
        %for j=1:N
            %G = [G11(i,j), G12(i,j); G12(i,j), G22(i,j)];
            %H(i,j) = det(G) - k*trace(G)^2;
            % Erweiterung um Rotationsinformation
            %winkel(i,j) = atan2(Iy_filt(i,j), Ix_filt(i,j));
        %end
    %end
    
    % Vectorized version of above code
    G_11 = reshape(G11,1,[]);
    G_22 = reshape(G22,1,[]);
    G_12 = reshape(G12,1,[]);
    G_det = G_11.*G_22 - G_12.^2;
    G_tr = G_11 + G_22;
    H = G_det - k*G_tr.^2;
    H = reshape(H,[M,N]);
    %Erweiterung um Winkel
    winkel = atan2(Iy_filt, Ix_filt);
    
    
    % Sinnvolle Behandlung der Merkmale, die n�her als ceil(segment_length/2) am Rand des Bildes sind
    
    % Eliminieren aller Merkmale, die kleiner als der Schwellwert tau sind
    corners = H;
    corners(corners < tau) = 0;
    %winkel(corners < tau) = 0;
    
    % finds the indices to the nonzero entries
    %[row,col] = find(corners);
    %merk = [col'; row'];
    
    %% Merkmalsvorbereitung (1.9)
    % Fuegen sie einen Nulland der Breite min_dist um die Matrix corners hinzu
    [M, N] = size(corners);
    c = zeros(M + 2*min_dist, N + 2*min_dist);
    c(min_dist+1:min_dist+M, min_dist+1:min_dist+N) = corners;
    
    % find nonzero entries
    [rowSub,colSub,v] = find(c);
    % Convert rowSub & colSub to linear matrix indices
    linearInd = sub2ind(size(c), rowSub, colSub);
    [~, I] = sort(v, 'descend');
    sorted_index = linearInd(I);
    
    %% Ueberschreibe corners variable mit c
    corners = c;    
    
    %% Harris Detector Kacheln (1.10)
    % Akkumulatorfeld
    % Legen Sie ein Akkumulatorfeld AKKA an (ein Eintrag pro Kachel)
    [rows,cols] = size(corners);
    tile_height = tile_size(1);
    tile_width = tile_size(2);
    AKKA = zeros(floor(rows/tile_height), floor(cols/tile_width));
    
    %[aRows,aCols] = size(AKKA);
    %m = zeros(2,size(sorted_index,1))    
    %m = zeros(2,N*aRows*aCols)
    %m = zeros(2,min(N*aRows*aCols, size(sorted_index,1)));
    %size(m);
    
    %% Merkmalsbestimmung mit Mindestabstand und Maximalzahl pro Kachel (1.11)
%     Cake = cake(min_dist);
%     tile_h = tile_size(1);
%     tile_w = tile_size(2);
%     %sizeAKKA = size(AKKA);    
%     %sizeSortedIndex = size(sorted_index)
%     merkmale = [];
%     
%     merk_index = 1;
%     % Finde die Kacheln zu sorted_index
%     for k=1:size(sorted_index,1)
%         % Berechne x,y Koordinaten dieses Merkmals
%         [y,x] = ind2sub(size(corners), sorted_index(k));
%         % Wenn dieses Feature noch nicht (durch Cake) auf 0 gesetzt wurde
%         if(corners(y,x) > 0)
%             % Berechne tile Indizes
%             j = ceil((x-min_dist) / tile_h);
%             i = ceil((y-min_dist) / tile_w);
%             if(AKKA(i,j) < N)
%                 % Erzeuge Cake um diese Koordinaten herum
%                 x_idx = (x-min_dist) : (x+min_dist);
%                 y_idx = (y-min_dist) : (y+min_dist);
%                 corners(y_idx, x_idx) = corners(y_idx, x_idx) .* Cake;
%                 % Fuege diesen Punkt den Merkmalen hinzu
%                 %merkmale(:, merk_index) = [x-min_dist; y-min_dist]; % x-min_dist; y-min_dist
%                 merkmale = [merkmale, [x-min_dist;y-min_dist]];
%                 merk_index = merk_index + 1;
%                 % Erh�he AKKA wert
%                 AKKA(i,j) = AKKA(i,j) + 1;
%             end
%         else
%             % Tue nichts
%         end
%     end
    
    % Get dims
    tile_height = tile_size(1);
    tile_width = tile_size(2);
    % Max merkmale = N per tile
    N_MAX = N*numel(AKKA);
    % Create linear index array
    corners_index = reshape(1:1:numel(corners), size(corners));
    cake_inv = ~cake(min_dist);
    % Reset merkmale
    merkmale = [];
    k = 1;
    
    % At most N_MAX features can be found
    while k < N_MAX
        %Stop if no more indices available
        if numel(sorted_index) == 0
            break
        end
        
        % Compute coordinates of highest index
        [row, col] = ind2sub(size(corners),sorted_index(1));
        % Compute tile of index and increase count
        tile_row = ceil((row-min_dist)/tile_height);
        tile_col = ceil((col-min_dist)/tile_width);
        AKKA(tile_row, tile_col) = AKKA(tile_row, tile_col) + 1;
        % Add merkmal
        merkmale = [merkmale, [col;row]];
        k = k + 1;
        % Always apply cake to index matrix
        row_idx = (row-min_dist):(row+min_dist);
        col_idx = (col-min_dist):(col+min_dist);
        tile = corners_index(row_idx,col_idx).*cake_inv;
        sorted_index(ismember(sorted_index,tile)) = [];
        % If limit is reached delete all tile indices from sorted_index
        if AKKA(tile_row, tile_col) == N
            % Delete all indices from tile that already have N features
            row_idx = 1+min_dist+(tile_row-1)*tile_height : min((min_dist+tile_row*tile_height), size(corners,1));
            col_idx = 1+min_dist+(tile_col-1)*tile_width : min((min_dist+tile_col*tile_width), size(corners,2));
            tile = corners_index(row_idx,col_idx);
            sorted_index(ismember(sorted_index,tile)) = [];
        end
    end
        merkmale = merkmale - min_dist;

    
    %% Plot (1.6)
    if(do_plot == true)
        figure
        imshow(input_image)
        hold on
        plot(merkmale(1,:), merkmale(2,:),'s')
        hold off
    end
end
