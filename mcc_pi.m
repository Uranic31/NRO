function [tocke_krog, tocke_kvadrat, zunaj_kroga] = mcc_pi(st_tock)
    % Generiranje naključnih točk znotraj kvadrata
    x = rand(1, st_tock);
    y = rand(1, st_tock);
    
    
    % Preverjanje, ali so točke znotraj kroga
    krog_noter = x.^2 + y.^2 <= 1;
    krog_zunaj = x.^2 + y.^2 > 1;

    % Koordinate točk znotraj kroga
    krog_x = x(krog_noter);
    krog_y = y(krog_noter);
    
    kvadrat_x = x;
    kvadrat_y = y;

    zunaj_kroga_x = x(krog_zunaj);
    zunaj_kroga_y = y(krog_zunaj);
    
    % Vrnemo koordinate točk znotraj kroga in kvadrata
    tocke_krog = [krog_x; krog_y];
    tocke_kvadrat = [kvadrat_x; kvadrat_y];
    zunaj_kroga = [zunaj_kroga_x; zunaj_kroga_y];
end
