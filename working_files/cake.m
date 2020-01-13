function Cake = cake(min_dist)
    % Die Funktion cake erstellt eine "Kuchenmatrix", die eine kreisfoermige
    % Anordnung von Nullen beinhaltet und den Rest der Matrix mit Einsen
    % auffuellt. Damit koennen, ausgehend vom staerksten Merkmal, andere Punkte
    % unterdrueckt werden, die den Mindestabstand hierzu nicht einhalten. 
    Cake = zeros(min_dist*2+1, min_dist*2+1);
    for i=1:min_dist*2+1
        for j = 1:min_dist*2+1
            if sqrt((i-min_dist-1)^2+(j-min_dist-1)^2) > min_dist
                Cake(i,j) = 1;
            end
        end
    end
    Cake = logical(Cake);
end
