clc;
clear all;

filename = './primer1mreza.txt';

X = [];
Y = [];
celice = [];
vozlisca_robnih_pogojev = {};
tipi_robnih_pogojev = [];
vrednosti_robnih_pogojev = [];
vrednosti_prestopa_toplote = [];

fid = fopen(filename);
prva_vrstica = fgetl(fid);
n_vozlisc = str2double(split(prva_vrstica));
n_vozlisc = n_vozlisc(2);

for line = 1:n_vozlisc
    line = fgetl(fid);   
    line = strrep(line, ';', ' ');
    line = strrep(line, ',', ' ');
    line_split = split(line);
    x = str2double(line_split(2));
    y = str2double(line_split(3));
    
    X(end+1) = x;
    Y(end+1) = y;
end

prazna_vrstica = fgetl(fid);

cell_line = fgetl(fid);
cells_string = split(cell_line);
n_celice = str2double(cells_string(2));

for line = 1:n_celice
    line = fgetl(fid);
    line = strrep(line, ',', ' ');
    line = strrep(line, ';', ' ');
    line_split = split(line);
    node1_id = str2double(line_split(2))+1;
    node2_id = str2double(line_split(3))+1;
    node3_id = str2double(line_split(4))+1;
    node4_id = str2double(line_split(5))+1;
    cell = [node1_id, node2_id, node3_id, node4_id];
    celice = [celice; cell];
end

prazna_vrstica = fgetl(fid);

robni_pogoji_line = fgetl(fid);
line_split = split(robni_pogoji_line);
n_pogoji = str2double(line_split(end));

for n = 1:n_pogoji
    line = fgetl(fid);
    tip_pogoja_vrstica = split(line);
    tip_pogoja = tip_pogoja_vrstica(3);
   
    if tip_pogoja == "temperatura"
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
        vrstica = split(line);
        koeficient_prestopa = str2double(vrstica(3));
        vrednosti_robnih_pogojev(end+1) = temperatura_prestopa;
        vrednosti_prestopa_toplote(end+1) = koeficient_prestopa;
    end
    
    line = fgetl(fid);
    st_vozlisc_v_robnem_pogoju = str2double(line);
    vozlisca_v_robnem_pogoju = [];
    for vozl = 1:st_vozlisc_v_robnem_pogoju
        line = fgetl(fid);
        id_vozlisce_v_robnem_pogoju = str2double(line)+1; 
        vozlisca_v_robnem_pogoju = [vozlisca_v_robnem_pogoju; id_vozlisce_v_robnem_pogoju];
    end
    vozlisca_robnih_pogojev{end+1}= vozlisca_v_robnem_pogoju; 
    
    line = fgetl(fid);
end

deltaX = 0.0125;
deltaY =0.0125;
k=24;

sosednja_vozlisca = [];

for node_id = 1:n_vozlisc
    node_i_neighbours = [-1,-1,-1,-1];
    
    for nd = 1:n_celice
        trenutna_celica = celice(nd, :);
        vozlisce1 = trenutna_celica(1);
        vozlisce2 = trenutna_celica(2);
        vozlisce3 = trenutna_celica(3);
        vozlisce4 = trenutna_celica(4);
        
        if (node_id == vozlisce1 | node_id == vozlisce2 | ... 
            node_id == vozlisce3 | node_id == vozlisce4 )
            
            for vozl = 1:4
                sosednje_vozlisce = trenutna_celica(vozl);
                if sosednje_vozlisce ~= node_id
                    x_obravnavano_vozl = X(node_id);
                    y_obravnavano_vozl = Y(node_id);
                    x_sosed = X(sosednje_vozlisce);
                    y_sosed = Y(sosednje_vozlisce);
                    
                    if (x_obravnavano_vozl-x_sosed) < 1e-9 && ...
                            (x_obravnavano_vozl-x_sosed) > -1e-9
                       
                        if (y_obravnavano_vozl-y_sosed) > 0
                            pozicija = 2;
                        else
                            pozicija = 4;
                        end
                        
                    elseif (y_obravnavano_vozl-y_sosed) < 1e-9 && ...
                            (y_obravnavano_vozl-y_sosed) > -1e-9
                   
                        if (x_obravnavano_vozl-x_sosed) > 0
                            pozicija = 1;
                        else
                            pozicija = 3;
                        end
                        
                    else
                        pozicija = -1;
                    end
                   
                    if pozicija ~= -1
                        node_i_neighbours(pozicija) = sosednje_vozlisce;
                    end
                end
                
            end
        end
        
    end
    sosednja_vozlisca = [sosednja_vozlisca; node_i_neighbours];
end

A = [];
for r = 1:n_vozlisc
    cols = [];
    for c = 1:n_vozlisc
        cols = [cols; 0];
    end
    A(end+1, :) = cols;
end

b = [];
for r = 1:n_vozlisc
    b = [b; 0];
end
 
for node_id = 1:n_vozlisc
    sosedi = sosednja_vozlisca(node_id, :);
    levi_sosed = sosedi(1);
    spodnji_sosed = sosedi(2);
    desni_sosed = sosedi(3);
    zgornji_sosed = sosedi(4);
    
    if (levi_sosed ~= -1 && spodnji_sosed ~=-1 && ...
        desni_sosed ~=-1 && zgornji_sosed ~=-1)
        A(node_id, levi_sosed) = 1;
        A(node_id, spodnji_sosed) = 1;
        A(node_id, desni_sosed) = 1;
        A(node_id, zgornji_sosed) = 1;
        A(node_id, node_id) = -4;
    else
        
        for robni_pogoj_id = 1:n_pogoji
            vozlisca_v_trenutnem_rp = cell2mat(vozlisca_robnih_pogojev(robni_pogoj_id));
       
            for id_vozlisce_rp = 1:size(vozlisca_v_trenutnem_rp)
                vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp(id_vozlisce_rp);
                if (node_id == vozlisce_v_trenutnem_rp)
                    tip_robnega_pogoja = tipi_robnih_pogojev(robni_pogoj_id);
                    vrednost = vrednosti_robnih_pogojev(robni_pogoj_id);
                    vrednost_prestopa = vrednosti_prestopa_toplote(robni_pogoj_id);
                end
            end
        end

        if (tip_robnega_pogoja==0)
            A(node_id, node_id) = 1;
            b(node_id) = vrednost;
        elseif (tip_robnega_pogoja==1)
            stevilo_sosedov = 0;
            for st = 1:4
                if (sosedi(st) ~= -1)
                    stevilo_sosedov = stevilo_sosedov + 1;
                end
            end
            
            if stevilo_sosedov == 3
                if levi_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id)-4;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
                if desni_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id) -4;
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
                if spodnji_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id)-4;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+2;
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+1;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
                if zgornji_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id)-4;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+2;
                    A(node_id, levi_sosed) = A(node_id, levi_sosed)+1;
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+1;
                    b(node_id) = -2*(vrednost*deltaX/k);
                end
            end
        elseif (tip_robnega_pogoja==2)
            stevilo_sosedov = 0;
            for st = 1:4
                if (sosedi(st) ~= -1)
                    stevilo_sosedov = stevilo_sosedov + 1;
                end
            end
            
            if stevilo_sosedov == 3
                if levi_sosed == -1
                    A(node_id, node_id) = A(node_id, node_id)-2*(vrednost_prestopa*deltaX/k+2);
                    A(node_id, desni_sosed) = A(node_id, desni_sosed)+2;
                    A(node_id, spodnji_sosed) = A(node_id, spodnji_sosed)+1;
                    A(node_id, zgornji_sosed) = A(node_id, zgornji_sosed)+1;
                    b(node_id) = b(node_id)-2*vrednost_prestopa*deltaX*vrednost/k;
                end
                if desni_sosed == -1
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

tic;
T=linsolve(A,b);

cas_resevanja_sistema = toc;




















%T = [];
%for T_zacetna = 1:n_vozlisc
 %   T(end+1) = 100;
%end

%for iitt = 1:1000
  %  for jj = 1:n_vozlisc
       % d = b(jj);

      %  for ii = 1:n_vozlisc
           % if(jj ~= ii)
               % d = d- A(jj, ii) * T(ii);
        
           % end
           % T(jj) = d / A(jj, jj);
       % end
    %end
%end

fileID = fopen('rezultat_vtk_linsolve.vtk','w');
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










