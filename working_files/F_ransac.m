function [Korrespondenzen_robust] = F_ransac(Korrespondenzen, varargin)
    % Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
    % robusten Korrespondenzpunktpaaren
    
    validProp = @(x) isnumeric(x) && (x>0) && (x<1);
    parser = inputParser;
    addOptional(parser,'epsilon',0.5,validProp);
    addOptional(parser,'p',0.5,validProp);
    addOptional(parser,'tolerance',0.01,@isnumeric);
    parse(parser,varargin{:});
    
    epsilon = parser.Results.epsilon;
    p = parser.Results.p;
    tolerance = parser.Results.tolerance;
    
    n = size(Korrespondenzen,2);
    

    x1_pixel = [ Korrespondenzen(1:2,:) ; ones(1,n) ];
    x2_pixel = [ Korrespondenzen(3:4,:) ; ones(1,n) ];


    %% RANSAC Algorithmus Vorbereitung
    k=8;
    s = log(1-p)/log(1-(1-epsilon)^k);
    largest_set_size = 0;
    largest_set_dist = inf;
    largest_set_F = zeros(3,3);


    %% RANSAC Algorithmus
    siz_korr = size(Korrespondenzen,2);
    for i=1:s
        idxp = randperm(siz_korr,k);
        F = achtpunktalgorithmus(Korrespondenzen(:,idxp));
        sd = sampson_dist(F, [Korrespondenzen(1:2,:);ones(1,siz_korr)], [Korrespondenzen(3:4,:);ones(1,siz_korr)]);
        idx_cons = find(sd < tolerance);
        n = size(idx_cons,2);
        dist = sum(sd(idx_cons),2);
        
        if n > largest_set_size
            largest_set_size = n;
            largest_set_F = F;
            largest_set_dist = dist;
            largest_set_idx = idx_cons;
        elseif ((n == largest_set_size) && (dist< largest_set_dist))
            largest_set_size = n;
            largest_set_F = F;
            largest_set_dist = dist;
            largest_set_idx = idx_cons;            
        end
        
    end
    Korrespondenzen_robust = Korrespondenzen(:,largest_set_idx);
    
    %Korrespondenzen_robust = {Korrespondenzen_robust, largest_set_F}
end