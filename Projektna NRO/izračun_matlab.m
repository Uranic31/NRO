clc;
clear all;

% Definiramo pot do datoteke z mrežo
filename = './primer1mreza.txt';


% Rešili bomo enostaven primer, ki je podan v obrazložitvi projektne
% naloge. Sicer skripta deluje za vse 4 primere, ki jih imate za projekt.
% Nastaviti morate pot v spremenljivki filename in pognati.
%
% Paziti moramo na indekse, saj se v C++ indeksiranje začne pri 0, v
% Matlabu pa pri 1. Tako bi mreža izgledala v Matlabu

% 1------2------3------4
% |      |      |      |
% |  1   |  2   |   3  |
% |      |      |      |
% 5------6------7------8
% |      |      |      |
% |  4   |  5   |   6  |
% |      |      |      |
% 9------10----11------12
% |      |      |      |
% |   7  |   8  |   9  |
% |      |      |      |
% 13-----14----15------16

% Imamo torej 16 vozlišč in 9 celic (včasih boste slišali tudi izraz 
% element=celica). Vsako vozlišče ima svoji koordinati x, y in svojo
% identifikacijsko številko. Vsaka celica ima svojo identifikacijsko
% številko, definirana pa je preko 4ih identifikacijskih številk vozlišč,
% ki ji pripadajo. Npr. za celico 1 imamo podano vozlišča 1,5,6,2, itd.
% Hkrati imamo podane robne pogoje:
% pogoj 1: T=50, vozlišča -> 1, 5, 9, 13
% pogoj 2: T=100, vozlišča -> 2, 3, 4
% pogoj 3: T=300, vozlišča -> 14, 15, 16
% pogoj 4: T_ext=200, h=1000, vozlišča -> 8, 12
%
% Vse te informacije si moramo nekako iz datoteke shraniti v program. Kako
% si bomo te strukture zamislili, je poljubno. Vse informacije lahko
% shranimo npr. v array-e (oz. v C++ v std::vector<>). Za vozlišča zapišemo
% vektorja X, Y, v katera shranimo koordinati. Za celice lahko uporabimo v
% c++ 'vector<vector<int>>'. Struktura v programu bi na koncu izgledala tako
%
% X = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
% Y = [3,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0]
% celice = [[1,5,6,2],
%           [2,6,7,3],
%           [3,7,8,4],
%           [5,9,10,6],
%           [6,10,11,7],
%           [7,11,12,8],
%           [9,13,14,10],
%           [10,14,15,11],
%           [11,15,16,12]]
%
% Do posameznih informacij dostopamo direktno preko indeksiranja. Npr. če
% želimo dobiti Y koordinato vozlišča 5, lahko rečemo
%
% y_vozlisca5 = Y[5] 
%
% Naslednji korak je shranjevanje robnih pogojev. Tu lahko npr. za svoj
% primer pogledate koliko robnih pogojev imate in jih direktno preberete,
% npr. če imate 4 robne pogoje, definirate spremenljivke
% robni_pogoj1_vozlisca, robni_pogoj1_temperatura in tako za vse ostale in
% shranite informacije direktno v te spremenljivke. V tem primeru bomo
% avtomatizirali proces tako, da bomo lahko prebrali katerokoli mrežo s
% poljubnimi robnimi pogoji. Lahko definiramo 4 vektorje
%
% vozlisca_robnih_pogojev;
% tipi_robnih_pogojev;
% vrednosti_robnih_pogojev;
% vrednosti_prestopa_toplote;
%
% V spremenljivko 'vozlisca_robnih_pogojev' bomo shranili vsa vozlišča,
% ki pripadajo posameznemu robnemu pogoju. 
%
% V spremenljivko 'tipi_robnih_pogojev' bomo shranili tip robnega pogoja, 
% ki je lahko v našem problemu: temperatura, toplotni tok ali prestop 
% toplote. Odločimo se, da bomo temperaturo označili z 0, toplotni tok z 1 
% in prestop toplote z 2.
%
% V spremenljivko vrednosti_robnih_pogojev bomo shranili temperaturo (robni
% pogoj tipa 0), toplotni tok (robni pogoj tipa 1) ali temperaturo okoliške
% tekočine T_ext (robni pogoj tipa 2).
%
% V spremenljivko 'vrednosti_prestopa_toplote' lahko shranimo koeficient
% prestopa toplote. Ta pride v poštev samo v primeru robnega pogoja 2. V
% ostalih dveh robnih pogojih, lahko sem dodamo NaN ali npr. -1.
%
% Za zgornji primer bi tako dobili
%
% vozlisca_robnih_pogojev = [[1, 5, 9, 13],
%                            [2, 3, 4],
%                            [14, 15, 16],
%                            [8, 12]]
%
% tipi_robnih_pogojev = [0, 0, 0, 2]
% vrednosti_robnih_pogojev = [50, 100, 300, 200]
% vrednosti_prestopa_toplote = [-1, -1, -1, 1000]
%
% Do posameznih informacij ponovno dostopamo preko
% indeksiranja. Npr. če želimo izvedeti tip za robni pogoj 3:
%
% tip_robni_pogoj3 = tipi_robnih_pogojev[3]
%
% S temi podatki bomo nato sestavili matriko za izračun.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prvi cilj je prebrati datoteko in iz nje izluščiti vse pomembne podatke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% potrebovali bomo torej sledeče parametre:
X = [];
Y = [];
celice = [];
vozlisca_robnih_pogojev = {};
tipi_robnih_pogojev = [];
vrednosti_robnih_pogojev = [];
vrednosti_prestopa_toplote = [];

% Sedaj lahko začnemo z branjem datoteke. Beremo vrstico po vrstico.
fid = fopen(filename);
% prva vrstica ima vedno podano število vseh vozlišč. S 'fgetl' preberemo
% naslednjo vrstico v datoteki. V C++ je identičen ukaz 'std::getline()'.
prva_vrstica = fgetl(fid);
n_vozlisc = str2double(split(prva_vrstica));
n_vozlisc = n_vozlisc(2);
% Za branje datoteke v C++ si poglejte tudi primer iz vaje 7, primer 2.
% Sicer lahko rečemo tudi
%
% std::ifstream file;
% file.open(filename);
% std::string str_tocke;
% indata >> str_tocke;
% std::string str_points;
% indata >> str_points;
% int n_nodes = std::stoi(str_points); // dobimo stevilo vozlisc
%
% Spremenljivka n_vozlisc nam pove, da imamo n vozlišč. To pomeni, da
% moramo prebrati naslednjih n vrstic, da dobimo koordinate X, Y.

for line = 1:n_vozlisc
    line = fgetl(fid);   
    line = strrep(line, ';', ' ');
    line = strrep(line, ',', ' ');
    line_split = split(line);
    x = str2double(line_split(2));
    y = str2double(line_split(3));
    % V tem primeru se najprej znebimo posebnih znakov ('special
    % character') ';' in ',' ter nato uporabimo split. Npr. v C++ lahko
    % rečemo:
    %
    % std::getline(file, s); // preberemo vrstico
    % std::replace(s.begin(), s.end(), ';', ' ');
    % std::replace(s.begin(), s.end(), ',', ' ');
    % std::istringstream iss(s);
    % int node_id;
    % double x;
    % double y;
    % iss >> node_id >> x>> y;
    
    X(end+1) = x;
    Y(end+1) = y;
    % V C++ uporabimo 'vector.push_back()'
end

% Ko končamo z branjem vseh koordinat, sledi prazna vrstica (poglejte
% datoteko z mrežo)
prazna_vrstica = fgetl(fid);

% Naslednja vrstica vsebuje podatke o celicah. Ponovimo isti princip kot
% pri branju točk
cell_line = fgetl(fid);
cells_string = split(cell_line);
n_celice = str2double(cells_string(2));

for line = 1:n_celice
    line = fgetl(fid);
    line = strrep(line, ',', ' ');
    line = strrep(line, ';', ' ');
    line_split = split(line);
    % Zaradi indeksiranja v Matlabu tukaj prištejem 1, v C++ ni potrebno
    node1_id = str2double(line_split(2))+1;
    node2_id = str2double(line_split(3))+1;
    node3_id = str2double(line_split(4))+1;
    node4_id = str2double(line_split(5))+1;
    cell = [node1_id, node2_id, node3_id, node4_id];
    celice = [celice; cell]; % V C++ uporabimo 'vector.push_back'
end

% Ko končamo s celicami, ponovno sledi prazna vrstica
prazna_vrstica = fgetl(fid);

% Nato sledijo robni pogoji
robni_pogoji_line = fgetl(fid);
line_split = split(robni_pogoji_line);
n_pogoji = str2double(line_split(end));
% V spremenljivki n_pogoji imamo število robnih pogojev v našem problemu.
% Sedaj bomo zapolnili vektorje, ki se nanašajo na robne pogoje
%vozlisca_robnih_pogojev = {};
%tipi_robnih_pogojev = [];
%vrednosti_robnih_pogojev = [];
%vrednosti_prestopa_toplote = [];

% Naredili bomo iteracijo čez vse robne pogoje in v vsaki iteraciji
% prebrali informacije za posamezni robni pogoj in jih shranili v naše
% vektorje. Tukaj lahko preberete za vsak robni pogoj posebej za vaš primer
% in ni potrebno posplošiti kot tukaj
for n = 1:n_pogoji
    line = fgetl(fid);
    tip_pogoja_vrstica = split(line);
    tip_pogoja = tip_pogoja_vrstica(3);
    % Ko dobimo tip pogoja, lahko uporabimo 'if' stavke, da preberemo
    % parametre in jih shranimo (tip robnega pogoja, vrednost temperature 
    % in vrednost prestopa)
    if tip_pogoja == "temperatura"
        % Če imamo temperaturo, preberemo naslednjo vrstico z vrednostjo in
        % shranimo v naše vektorje
        tipi_robnih_pogojev(end+1) = 0;
        line = fgetl(fid);
        vrstica = split(line);
        temperatura = str2double(vrstica(2));
        vrednosti_robnih_pogojev(end+1) = temperatura;
        vrednosti_prestopa_toplote(end+1) = -1;
    elseif tip_pogoja == "toplotni"
        tipi_robnih_pogojev(end+1) = 1;
        line = fgetl(fid);
        vrstica = split(line);
        toplotni_tok = str2double(vrstica(3));
        vrednosti_robnih_pogojev(end+1) = toplotni_tok;
        vrednosti_prestopa_toplote(end+1) = -1;
    elseif tip_pogoja == "prestop"
        tipi_robnih_pogojev(end+1) = 2;
        line = fgetl(fid);
        vrstica = split(line);
        temperatura_prestopa = str2double(vrstica(2));
        line = fgetl(fid);
        % pri prestopu moramo prebrati še eno vrstico, saj imamo dve
        % informaciji, poleg okoliške temperature še prestop
        vrstica = split(line);
        koeficient_prestopa = str2double(vrstica(3));
        vrednosti_robnih_pogojev(end+1) = temperatura_prestopa;
        vrednosti_prestopa_toplote(end+1) = koeficient_prestopa;
    end
    
    % shranimo še ID-je vozlišč, ki pripadajo robnemu pogoju
    line = fgetl(fid);
    st_vozlisc_v_robnem_pogoju = str2double(line);
    vozlisca_v_robnem_pogoju = [];
    for vozl = 1:st_vozlisc_v_robnem_pogoju
        line = fgetl(fid);
        % Prištejemo 1 zaradi indeksiranja v Matlabu
        id_vozlisce_v_robnem_pogoju = str2double(line)+1; 
        vozlisca_v_robnem_pogoju = [vozlisca_v_robnem_pogoju; id_vozlisce_v_robnem_pogoju];
    end
    vozlisca_robnih_pogojev{end+1}= vozlisca_v_robnem_pogoju; 
    
    % Na koncu preberemo tudi naslednjo vrstico, ki je prazna. Tako v novi
    % iteraciji takoj skočimo v vrstico, kjer imamo definiran robni pogoj
    line = fgetl(fid);
end

% Sedaj smo prebrali vse informacije, ki jih potrebujemo. Na tem mestu naše
% spremenljivke izgledajo takole

% X = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
% Y = [3,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0]
% celice = [[1,5,6,2],
%           [2,6,7,3],
%           [3,7,8,4],
%           [5,9,10,6],
%           [6,10,11,7],
%           [7,11,12,8],
%           [9,13,14,10],
%           [10,14,15,11],
%           [11,15,16,12]]
% vozlisca_robnih_pogojev = {[1, 5, 9, 13],
%                            [2, 3, 4],
%                            [14, 15, 16],
%                            [8, 12]}
%
% tipi_robnih_pogojev = [0, 0, 0, 2]
% vrednosti_robnih_pogojev = [50, 100, 300, 200]
% vrednosti_prestopa_toplote = [-1, -1, -1, 1000]

% 2.del
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Priprava matrike A in vektorja b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iskanje sosedov vozlišč
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tukaj za vaš primer definirajte koeficient prevoda toplote. Tudi za vaš
% primer je deltaX in deltaY konstantna. Lahko jih definirate že na začetku
% skripte
deltaX = 1;
deltaY = 1;
k=1;

% V vektor sosednja vozlišča bomo za vsako vozlišče shranili
% identifikacijske točke sosednjih vozlišč. Potrebujemo vozlišča, ki se
% nahajajo levo, desno, zgoraj in spodaj od obravnavanega vozlišča.
sosednja_vozlisca = [];
% npr. v C++ bi definirali:
%
% vector<vector<int>> sosednja_vozlisca;
%
% ter nato dodajamo z uporabo funkcije 'push_back()'.

% za nastavitev formule končnih razlik za vsako vozlišče torej potrebujemo
% sosednja vozlišča, ki se nahajajo levo, desno, zgoraj in spodaj. 
% V ta namen lahko naredimo 'for' zanko skozi vsa vozlišča in napišemo 
% algoritem, ki nam poda omenjena sosednja vozlišča. Težava se bo pojavila
% pri vozliščih, ki se nahajajo na robu in nimajo vseh štirih sosedov.
% Ponovno lahko inicializiramo vektor
%
% node_i_neighbours = [-1,-1,-1,-1];
%
% in predpostavimo interno pravilo:
% levi sosed -> prvi element v vektorju
% spodnji sosed -> drugi element v vektorju
% desni sosed -> tretji element v vektorju
% zgornji sosed -> četrti element v vektorju
% 
% Če je vrednost -1, potem na tistem položaju ni soseda
% 
% poglejmo si nekaj primerov za naš primer:
% 1------2------3------4
% |      |      |      |
% |  1   |  2   |   3  |
% |      |      |      |
% 5------6------7------8
% |      |      |      |
% |  4   |  5   |   6  |
% |      |      |      |
% 9------10----11------12
% |      |      |      |
% |   7  |   8  |   9  |
% |      |      |      |
% 13-----14----15------16
%
% vozlišče 6:
% node_i_neighbours = [5,10,7,2];
%
% vozlišče 14:
% node_i_neighbours = [13,-1,15,10];
%
% vozlišče 16:
% node_i_neighbours = [15,-1,-1,12];

% Naredimo iteracijo po vozliščih
for node_id = 1:n_vozlisc
    % Spremenljivka node_id je identifikacijska številka trenutno 
    % obravnavanega vozlišča
    % V spremenljivko node_i_neighbours bomo shranili sosede trenutno obravnavanega
    % vozlišča. Na koncu bomo to spremenljivko dodali v vektor 
    % sosednja_vozlisca, da bomo imeli vse shranjeno.
    node_i_neighbours = [-1,-1,-1,-1];
    
    % Za posamezno vozlišče naredimo iteracijo po vseh celicah in preverimo
    % ali se obravnavano vozlišče nahaja v celici. Če se, potem gremo
    % preverjati, kaj je levo, desno, ...
    for nd = 1:n_celice
        trenutna_celica = celice(nd, :);
        % V C++ lahko rečete
        %
        % vector<int> trenutna_celica = celice[nd];
        %
        vozlisce1 = trenutna_celica(1);
        vozlisce2 = trenutna_celica(2);
        vozlisce3 = trenutna_celica(3);
        vozlisce4 = trenutna_celica(4);
        
        % sedaj moramo ugotoviti ali se obravnavano vozlišče nahaja
        % v celici
        if (node_id == vozlisce1 | node_id == vozlisce2 | ... 
            node_id == vozlisce3 | node_id == vozlisce4 )
            % če se nahaja v celici, pomeni, da se od treh preostalih
            % vozlišč dve nahajata levo, desno, zgoraj ali spodaj. Ti dve
            % vozlišči moramo identificirati, zato da lahko nastavimo
            % formulo za končne razlike
            
            for vozl = 1:4
                % naredimo iteracijo po trenutni celici in preverimo ali
                % je ID enak obravnavanemu vozlišču
                sosednje_vozlisce = trenutna_celica(vozl);
                if sosednje_vozlisce ~= node_id
                    % v kolikor ni enak, moramo preveriti, ali se sosednje 
                    % vozlišče glede na trenutno obravnavano 
                    % nahaja levo, desno, spodaj, zgoraj ali na diagonali
                    
                    % To lahko preverimo npr. tako da pogledamo ali
                    % sta vozlišči poravnani vertikalno ali horizontalno
                    
                    % izpišemo koordinati obeh vozlišč
                    x_obravnavano_vozl = X(node_id);
                    y_obravnavano_vozl = Y(node_id);
                    x_sosed = X(sosednje_vozlisce);
                    y_sosed = Y(sosednje_vozlisce);
                    
                    % priporočljivo je, da za to preverjanje napišete
                    % funkcijo in jo pokličete, da je koda bolj pregledna
                    % npr. C++:
                    %
                    % int preveri(double x_trenutno, double y_trenutno,
                    %             x_sosed, y_sosed){...}
                    %
                    
                    if (x_obravnavano_vozl-x_sosed) < 1e-9 && ...
                            (x_obravnavano_vozl-x_sosed) > -1e-9
                        % sosednje vozlišče se nahaja vertikalno
                        % preverimo ali spodaj - zgoraj
                        if (y_obravnavano_vozl-y_sosed) > 0
                            % se nahaja spodaj - 2
                            pozicija = 2;
                        else
                            % se nahaja zgoraj - 4
                            pozicija = 4;
                        end
                        
                    elseif (y_obravnavano_vozl-y_sosed) < 1e-9 && ...
                            (y_obravnavano_vozl-y_sosed) > -1e-9
                        % sosednje vozlišče se nahaja horizontalno
                        % preverimo ali levo - desno
                        if (x_obravnavano_vozl-x_sosed) > 0
                            % se nahaja levo - 1
                            pozicija = 1;
                        else
                            % se nahaja desno - 3
                            pozicija = 3;
                        end
                        
                    else
                        % Ce vozlišči nista poravnani vertikalno ali
                        % horizontalno ne storimo ničesar
                        pozicija = -1;
                    end
                   
                    % v kolikor smo dobili eno izmed sosednjih vrednosti
                    % 1, 2, 3, 4
                    if pozicija ~= -1
                        node_i_neighbours(pozicija) = sosednje_vozlisce;
                    end
                end
                
            end
        end
        
    end
    % Na koncu dodamo vektor sosednjih vozlišč, v 'sosednja_vozlisca', da
    % imamo nabor sosedov za vsako vozlišče
    sosednja_vozlisca = [sosednja_vozlisca; node_i_neighbours];
end

% Sedaj smo zaključili z iskanjem sosedov. Informacije, ki jih bomo
% potrebovali, so sedaj zbrane v vektorju sosednja_vozlisca, ki izgleda 
% takole
%
% sosednja_vozlisca = [[-1,     5,     2,    -1],
%		               [ 1,     6,     3,    -1],
%		               [ 2,     7,     4,    -1],
%		               [ 3,     8,    -1,    -1],
%		               [-1,     9,     6,     1],
%		               [ 5,    10,     7,     2],
%		               [ 6,    11,     8,     3],
%		               [ 7,    12,    -1,     4],
%		               [-1,    13,    10,     5],
%		               [ 9,    14,    11,     6],
%		               [10,    15,    12,     7],
%		               [11,    16,    -1,     8],
%		               [-1,    -1,    14,     9],
%		               [13,    -1,    15,    10],
%		               [14,    -1,    16,    11],
%		               [15,    -1,    -1,    12]]
%
% Do izpisa sosednjih vozlišč za posamezno vozlišče dostopamo kar preko
% indeksiranja. Npr. če želimo seznamo sosednjih vozlišč za vozlišče 5
%
% sosedi_za_vozl_5 = sosednja_vozlisca(5, :) 
%
% Oz. v C++
%
% vector<int> sosedi_za_vozl_5 = sosednja_vozlisca[5];
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradnja matrike A in vektorja b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Najprej inicializiramo matriko A z vrednostmi 0. V Matlabu bi lahko 
% uporabili built-in funkcijo zeros(). V C++ lahko naredimo iteracijo. 
% Naredimo jo tudi v tem primeru za Matlab. Matrika A ima toliko vrstic kot
% je vozlišč. 
A = [];
for r = 1:n_vozlisc
    cols = [];
    for c = 1:n_vozlisc
        cols = [cols; 0];
    end
    A(end+1, :) = cols;
end

% Nato inicializiramo vektor b, ki ima toliko elementov kot je vozlišč
b = [];
for r = 1:n_vozlisc
    b = [b; 0];
end

% Sedaj ko imamo za vsako vozlišče definirane tudi njegove sosede, lahko 
% ponovno naredimo iteracijo skozi vozlišča in zapišemo člene v matriko. S
% pomočjo matrike sosednja_vozlisca kjer imamo podana sosednja vozlišča za 
% vsako vozlišče

% Naredimo iteracijo skozi vrstice 
for node_id = 1:n_vozlisc
    % node_id je ID obravnavanega vozlišča in je enak ID-ju vrstice v A
    % izpišemo vse sosede, ki jih ima obravnavano vozlišče
    sosedi = sosednja_vozlisca(node_id, :);
    levi_sosed = sosedi(1);
    spodnji_sosed = sosedi(2);
    desni_sosed = sosedi(3);
    zgornji_sosed = sosedi(4);
    
    % Sedaj je potrebno ugotoviti, ali je vozlišče na robu in če je, kakšen
    % je njegov robni pogoj
    % V kolikor nobena vrednost ni -1, vozlišče ni na robu domene
    if (levi_sosed ~= -1 && spodnji_sosed ~=-1 && ...
        desni_sosed ~=-1 && zgornji_sosed ~=-1)
        % nismo na robu, zato zapisemo enačbo (12) iz predloge za projekt
        A(node_id, levi_sosed) = 1;
        A(node_id, spodnji_sosed) = 1;
        A(node_id, desni_sosed) = 1;
        A(node_id, zgornji_sosed) = 1;
        A(node_id, node_id) = -4;
    else
        % če ima eno izmed sosednjih vozlišč -1, pomeni, da je obravnavano
        % vozlišče na robu, kar pomeni, da ima predpisan robni pogoj.
        % Najprej moramo ugotoviti, kateremu robnemu pogoju pripada
        % obravnavano vozlišče, nato pa sestavimo enačbo in jo damo v
        % matriko A in vektor b
        
        % Za to, kateremu robnemu pogoju pripada vozlišče, bi bilo dobro
        % napisati funkcijo, da je koda bolj pregledna. Tukaj bomo
        % nadaljevali brez funkcije. Napisali bomo iteracijo skozi vse
        % robne pogoje ter nato skozi vsa vozlišča posameznega robnega
        % pogoja in preverili, ali se obravnavano vozlišče tam nahaja
        
        for robni_pogoj_id = 1:n_pogoji
            vozlisca_v_trenutnem_rp = cell2mat(vozlisca_robnih_pogojev(robni_pogoj_id));
            % Matlab ne pusti definicije vektorja, ki bi imel shranjene
            % različne dolžine vektorjev, zato uporabimo tip strukture cell
            % in pretvorimo s cell2mat. V C++ to ni težava in izgleda
            % takole:
            %
            % vector<int> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[robni_pogoj_id];
            %        
            % naredimo iteracijo skozi vsa vozlišča v trenutnem robnem
            % pogoju in pogledamo, če je obravnavano vozlišče tam
            for id_vozlisce_rp = 1:size(vozlisca_v_trenutnem_rp)
                vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp(id_vozlisce_rp);
                if (node_id == vozlisce_v_trenutnem_rp)
                    % če je obravnavano vozlišče v skupini vozlišč, ki 
                    % pripadajo robnemu pogoju, izpišemo parametre za ta
                    % robni pogoj
                    
                    tip_robnega_pogoja = tipi_robnih_pogojev(robni_pogoj_id);
                    % tip robnega pogoja je lahko 0 (temperatura), 1
                    % (toplotni tok) ali 2 (prestop toplote)
                    
                    vrednost = vrednosti_robnih_pogojev(robni_pogoj_id);
                    % vrednost je temperatura (za tip 0), toplotni tok (za
                    % tip 1) ali okoliska temperatura (za tip 2)
                    
                    vrednost_prestopa = vrednosti_prestopa_toplote(robni_pogoj_id);
                    % vrednost prestopa je -1 (za tip 0), -1 (za
                    % tip 1) ali koeficient prestopa (za tip 2). Uporabimo
                    % samo pri tipu 2
                end
            end
        end
        
        % Sedaj ko smo šli čez vse robne pogoje in ugotovili, katerega
        % potrebujemo, lahko zapisemo enačbe glede na to katere sosede ima
        % obravnavano vozlišče na robu
        
        if (tip_robnega_pogoja==0)
            % imamo temperaturo, ki je določena z vrednostjo.
            % vrednost_prestopa je -1 in je tu ne potrebujemo. Vseeno je,
            % ali je vozlišče brez dveh sosedov 
            % (npr. v kotu -> nima npr. levega in zgornjega soseda) 
            % ali pa je vozlišče brez enega soseda (nima npr. levega soseda)
            A(node_id, node_id) = 1;
            b(node_id) = vrednost;
        elseif (tip_robnega_pogoja==1)
            % imamo toplotni tok. Tu je enačba odvisna od tega, koliko
            % sosedov ima vozlišče in kateri so ti sosedi (levo, desno, ...)
            % Zato najprej izračunamo, koliko sosedov ima vozlišče
            stevilo_sosedov = 0;
            for st = 1:4
                if (sosedi(st) ~= -1)
                    stevilo_sosedov = stevilo_sosedov + 1;
                end
            end
            
            % Sedaj lahko na podlagi števila sosedov nastavimo enačbe
            if stevilo_sosedov == 3
                if levi_sosed == -1
                    % Če levo od obravnavanega vozlišča ni soseda, pomeni,
                    % da je vozlišče na levem robu, uporabimo
                    A(node_id, node_id) = A(node_id, node_id)-4;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
                if desni_sosed == -1
                    % Če desno od obravnavanega vozlišča ni soseda, pomeni,
                    % da je vozlišče na desnem robu, uporabimo enačbo (11)
                    % iz predloge za projekt
                    A(node_id, node_id) = A(node_id, node_id) -4;
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
                if spodnji_sosed == -1
                    % če spodaj od obravnavanega vozlišča ni soseda, pomeni,
                    % da je vozlišče na spodnjem robu
                    A(node_id, node_id) = A(node_id, node_id)-4;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+2;
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+1;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
                if zgornji_sosed == -1
                    % če zgoraj od obravnavanega vozlišča ni soseda, pomeni,
                    % da je vozlišče na zgornjem robu
                    A(node_id, node_id) = A(node_id, node_id)-4;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+2;
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+1;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
            end
        elseif (tip_robnega_pogoja==2)
            % imamo prestop toplote. Tu je enačba ponovno odvisna od tega, koliko
            % sosedov ima vozlišče in kateri so ti sosedi (levo, desno, ...)
            % Zato najprej izračunamo, koliko sosedov ima vozlišče
            stevilo_sosedov = 0;
            for st = 1:4
                if (sosedi(st) ~= -1)
                    stevilo_sosedov = stevilo_sosedov + 1;
                end
            end
            
            % Sedaj lahko na podlagi števila sosedov nastavimo enačbe. Ni
            % potrebno nastaviti za vse kombinacije, lahko samo za vaš
            % primer
            if stevilo_sosedov == 3
                % lahko tudi preverimo, kateri sosed manjka in zapišemo
                % enačbo
                if levi_sosed == -1
                    % če levo od obravnavanega vozlišča ni soseda, pomeni,
                    % da je vozlišče na levem robu
                    A(node_id, node_id) = A(node_id, node_id)-2*(vrednost_prestopa*deltaX/k+2);
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = b(node_id)-2*vrednost_prestopa*deltaX*vrednost/k;
                end
                if desni_sosed == -1
                    % uporabimo enačbo 9 iz predloge
                    A(node_id, node_id) = A(node_id, node_id)-2*(vrednost_prestopa*deltaX/k+2);
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = b(node_id)-2*vrednost_prestopa*deltaX*vrednost/k;
                end
                if spodnji_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id)-2*(vrednost_prestopa*deltaX/k+2);
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+1;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+2;
                    b(node_id) = -2*vrednost_prestopa*deltaX*vrednost/k;
                end
                if zgornji_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id)-2*(vrednost_prestopa*deltaX/k+2);
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+1;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+1;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+2;
                    b(node_id) = -2*vrednost_prestopa*deltaX*vrednost/k;
                end
            end
        end
    end
end

% Sedaj imamo sestavljeno matriko A in vektor b, ki sta enaka kot sta
% zapisani v predlogi za projekt.

% 3. korak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rešitev sistema enačb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sedaj smo sestavili matriko A in b. Naslednji korak je rešitev sistema
% enačb. V principu se v Matlabu uporabljajo built-in funkcije, kot so
% linsolve,...

% V C++ se lahko uporabijo funkcije iz knjižnice Eigen, kar bo pokazano na 
% vajah/predavanjih. Za relativno enostaven primer sistema enačb pri 
% projektu, kjer je večina elementov blizu diagonale in imamo redko matriko 
% (ang. 'sparse matrix') lahko uporabimo tudi Gauss-Seidel metodo za sistem 
% enačb

% Najprej inicializiramo vektor rešitev na začetno vrednost. Kot začetno
% vrednost lahko izberemo eno izmed temperatur pri robnem pogoju
T = [];
for T_zacetna = 1:n_vozlisc
    T(end+1) = 100;
end

% Nato izberemo večje število iteracij ali določimo red napake in v vsaki 
% iteraciji izračunamo vrednosti v T. Sledimo psevdokodi v predlogi za 
% projekt
for iitt = 1:1000
    % for zanka za izračun T
    for jj = 1:n_vozlisc
        d = b(jj);

        for ii = 1:n_vozlisc
            if(jj ~= ii)
                d = d- A(jj, ii) * T(ii);
        
            end
            T(jj) = d / A(jj, jj);
        end
    end
end

% Dobimo rešitev temperatur v obliki vektorja T.
%
% T = [50,100,100,100,50,118.74,
%      156.22,199.90,50,
%      168.75,206.25,
%      200.05,50,300,300,300]
%

% 4. korak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vizualizacija podatkov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Za vizualizacijo je treba podatke shraniti v primeren format in jih
% odpreti s primernim grafičnim vmesnikom. V Matlabu lahko uporabimo 
% built-in funkcije za risanje grafov, itd..

% Na prihodnjih vajah si bomo pogledali VTK in ParaView. V C++ obstajajo
% knjižnice, ki avtomatično shranijo mreže in rezultate v ustrezen VTK
% format in si jih bomo pogledali kasneje na vajah/predavanjih. Tukaj je
% primer ročnega zapisa shranjevanja v datoteko VTK.

fileID = fopen('rezultat_vtk.vtk','w');
fprintf(fileID,'# vtk DataFile Version 3.0\n');
fprintf(fileID,'Mesh_1\n');
fprintf(fileID,'ASCII\n');
fprintf(fileID,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fileID,'POINTS %d float\n', n_vozlisc);
for koordinata_id = 1:n_vozlisc
    fprintf(fileID,'%f %f 0\n', X(koordinata_id), Y(koordinata_id));
end    
fprintf(fileID, '\n');
fprintf(fileID,'CELLS %d %d\n', n_celice, n_celice*5);
for celica_id = 1:n_celice
    vozl_id1 = celice(celica_id, 1)-1;
    vozl_id2 = celice(celica_id, 2)-1;
    vozl_id3 = celice(celica_id, 3)-1;
    vozl_id4 = celice(celica_id, 4)-1;
    fprintf(fileID,'4 %d %d %d %d\n', vozl_id1, vozl_id2, vozl_id3, vozl_id4);
end   
fprintf(fileID, '\n');
fprintf(fileID,'CELL_TYPES %d\n', n_celice);
for celica_id = 1:n_celice
    fprintf(fileID, '9\n');
end
fprintf(fileID, '\n');
fprintf(fileID, 'POINT_DATA %d\n', n_vozlisc);
fprintf(fileID, 'SCALARS Temperature float 1\n', n_vozlisc);
fprintf(fileID, 'LOOKUP_TABLE default\n', n_vozlisc);
for koordinata_id = 1:n_vozlisc
    fprintf(fileID,'%f\n', T(koordinata_id));
end   










