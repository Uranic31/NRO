function calc_pi()
    st_tock=30000;
    R = 1;
    [v_krogu, v_kvadratu, zunaj]=mcc_pi(st_tock);
    [ocena_pi, napaka_pi]=area_pi(v_krogu, v_kvadratu);
    visualize(v_krogu, zunaj, R);

    fprintf('Ocenjena vrednost π: %f\n', ocena_pi);
    fprintf('Napaka: %f\n', napaka_pi);
end

function [ocena,napaka]=area_pi(krog,kvadrat)
    ocena = 4 * size(krog, 2) / size(kvadrat, 2);
    realni = pi;
    napaka = abs(ocena - realni);
end

function visualize(kro,kvad,r)
    figure;
    hold on;
    scatter(kro(1,:),kro(2,:),'red','x');
    scatter(kvad(1,:),kvad(2,:),'yellow','filled');
    izrisi_kroznico(r);
    hold off;
end

 function izrisi_kroznico(r)
    theta = linspace(0, 2 * pi, 1000);
    kroznica = @(r, theta) [r * cos(theta); r * sin(theta)];
    kroznica_tocke = kroznica(r, theta);
    plot(kroznica_tocke(1, :), kroznica_tocke(2, :), 'green',LineWidth=2);
    axis equal;
    axis([0, 1, 0, 1]);
 end

