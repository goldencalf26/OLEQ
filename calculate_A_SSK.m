function A_SSK = calculate_A_SSK(psi, alph, phi)
    % из ОСК в ССК
    Apsi = [1, 0, 0; 0 cos(psi), sin(psi); 0, -sin(psi), cos(psi)];
    Aalph = [cos(alph), 0, -sin(alph); 0, 1, 0; sin(alph), 0, cos(alph)];
    Apfi = [1, 0, 0; 0 cos(phi), sin(phi); 0, -sin(phi), cos(phi)];

    A_SSK = Apfi*Aalph*Apsi;
end