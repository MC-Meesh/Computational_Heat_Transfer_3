% Heat Transfer - Spring 2023 Computational Assignment 3
% Michael Allen, Cullen Hirstius, Hayden Payne
% 4/21/23

clc; clear; close all;

%% Scenario 1
%gap in between pot and stove
heatTransferArray1 = [];
viewFactorArray = [];
heightArray = [];

for H = .01:.01:.3
    heightArray = [heightArray, H];
    %problem setup
    
    %both pot and stove can be treated as opaque (no transmissivity), gray
    %(epsilon = alpha), and diffuse (black body radiation)
    
    diameter_pot = .2; %diameter of the pot [m]
    T_pot = 373; %temperature of pot [K]
    epsilon_pot = .3; %emmisivity of pot
    A_pot = pi * (diameter_pot / 2)^2; % area of bottom surface of pot [m^2]
    
    
    diameter_stove = diameter_pot; %diameter of stovetop [m]
    T_stove = 1273; %temperature of stovetop [K]
    epsilon_stove = 1; %emmisivity of stovetop, blackbody
    A_stove = pi * (diameter_stove / 2)^2; % area of stove of pot [m^2]
    
    T_surr = 27+273; %temperature of surroundings [K]
    
    A_surr = pi*diameter_pot * H;
    
    
    % Calculate View Factors
    R_j = (diameter_pot / 2)/H;
    R_i = (diameter_stove / 2)/H;
    S = 1 + (1+R_j^2) / R_i^2;
    F_ij = .5*(S-(S^2-4*(R_j/R_i)^2)^(.5));
    
    
    %stove view factors
    F_stove_pot = F_ij;
    F_stove_stove = 0;
    F_stove_surr = 1 - F_stove_pot;
    
    %pot view factors
    F_pot_pot = 0;
    F_pot_stove = (A_stove / A_pot) * F_stove_pot;
    F_pot_surr = 1 - F_pot_stove;
    
    %side view factors
    F_surr_pot = (A_pot / A_surr) * F_pot_surr;
    F_surr_stove = (A_stove / A_surr) * F_stove_surr;

    viewFactorArray = [viewFactorArray, F_stove_pot];
    
    
    % Setup System for Finding q_pot
    sympref('FloatingPointOutput', true);
    
    sigma = 5.67*10^(-8);
    
    syms J_1 J_2 J_3;
    stove1 = J_1 == sigma*T_stove^4;
    pot2 = (sigma*T_pot^4 - J_2)*(epsilon_pot * A_pot) / (1-epsilon_pot) == A_pot * (F_pot_stove*(J_2-J_1) + F_pot_surr*(J_2-J_3));
    surroundings3 = J_3 == (sigma)*T_surr^4;
    soln = solve([stove1 pot2 surroundings3], [J_1 J_2 J_3]);

    if H > .099 && H < .101
       disp("Scenario 1 Radioosity Values");
       disp([soln.J_1 soln.J_2 soln.J_3]);
    end
    
    q_pot = (sigma*T_pot^4 - soln.J_2)*(epsilon_pot * A_pot) / (1-epsilon_pot);
    heatTransferArray1 = [heatTransferArray1, q_pot];
end

%% Scenario 2
%assume Al foil acting as cylindrical shroud

heatTransferArray2 = [];

for H = .01:.01:.3
    %problem setup
    
    %both pot and stove can be treated as opaque (no transmissivity), gray
    %(epsilon = alpha), and diffuse (black body radiation)
    
    diameter_pot = .2; %diameter of the pot [m]
    T_pot = 373; %temperature of pot [K]
    epsilon_pot = .3; %emmisivity of pot
    A_pot = pi * (diameter_pot / 2)^2; % area of bottom surface of pot [m^2]
    
    
    diameter_stove = diameter_pot; %diameter of stovetop [m]
    T_stove = 1273; %temperature of stovetop [K]
    epsilon_stove = 1; %emmisivity of stovetop, blackbody
    A_stove = pi * (diameter_stove / 2)^2; % area of stove of pot [m^2]

    T_surr = 27+273; %temperature of surroundings [K]

    A_foil = pi*diameter_pot * H; %area of foil replaced that of virtual surface from problem 1
    T_foil = T_surr; %foil is thin, so T_foil is same of surroundings
    epsilon_foil = .1; %emissivity of Al foil
    
    
    
    
    
    % Calculate View Factors
    R_j = (diameter_pot / 2)/H;
    R_i = (diameter_stove / 2)/H;
    S = 1 + (1+R_j^2) / R_i^2;
    F_ij = .5*(S-(S^2-4*(R_j/R_i)^2)^(.5));
    
    
    %stove view factors
    F_stove_pot = F_ij;
    F_stove_stove = 0;
    F_stove_foil = 1 - F_stove_pot;
    
    %pot view factors
    F_pot_pot = 0;
    F_pot_stove = (A_stove / A_pot) * F_stove_pot;
    F_pot_foil = 1 - F_pot_stove;
    
    %side view factors
    F_foil_pot = (A_pot / A_foil) * F_pot_foil;
    F_foil_stove = (A_stove / A_foil) * F_stove_foil;
    
    
    % Setup System for Finding q_pot
    sympref('FloatingPointOutput', true);
    
    sigma = 5.67*10^(-8);
    
    syms J_1 J_2 J_3;
    stove1 = J_1 == sigma*T_stove^4;
    pot2 = (sigma*T_pot^4 - J_2)*(epsilon_pot * A_pot) / (1-epsilon_pot) == A_pot * (F_pot_stove*(J_2-J_1) + F_pot_foil*(J_2-J_3));
    foil3 = (sigma*T_foil^4 - J_3)*(epsilon_foil * A_foil) / (1-epsilon_foil) == A_foil * (F_foil_stove*(J_3-J_1) + F_foil_pot*(J_3-J_2));
    soln = solve([stove1 pot2 foil3], [J_1 J_2 J_3]);


     if H > .099 && H < .101
       disp("Scenario 2 Radioosity Values");
       disp([soln.J_1 soln.J_2 soln.J_3]);
    end
    
    q_pot = (sigma*T_pot^4 - soln.J_2)*(epsilon_pot * A_pot) / (1-epsilon_pot);
    heatTransferArray2 = [heatTransferArray2, q_pot];
end

%% Scenario 3
%Assume shroud is made of black paint

heatTransferArray3 = [];

for H = .01:.01:.3
    %problem setup
    
    %both pot and stove can be treated as opaque (no transmissivity), gray
    %(epsilon = alpha), and diffuse (black body radiation)
    
    diameter_pot = .2; %diameter of the pot [m]
    T_pot = 373; %temperature of pot [K]
    epsilon_pot = .3; %emmisivity of pot
    A_pot = pi * (diameter_pot / 2)^2; % area of bottom surface of pot [m^2]
    
    
    diameter_stove = diameter_pot; %diameter of stovetop [m]
    T_stove = 1273; %temperature of stovetop [K]
    epsilon_stove = 1; %emmisivity of stovetop, blackbody
    A_stove = pi * (diameter_stove / 2)^2; % area of stove of pot [m^2]

    T_surr = 27+273; %temperature of surroundings [K]

    A_shroud = pi*diameter_pot * H; %area of foil replaced that of virtual surface from problem 1
    T_shroud = T_surr; %foil is thin, so T_shroud is same of surroundings
    epsilon_shroud = .9; %emissivity of Al shroud
    
    
    
    
    
    % Calculate View Factors
    R_j = (diameter_pot / 2)/H;
    R_i = (diameter_stove / 2)/H;
    S = 1 + (1+R_j^2) / R_i^2;
    F_ij = .5*(S-(S^2-4*(R_j/R_i)^2)^(.5));
    
    
    %stove view factors
    F_stove_pot = F_ij;
    F_stove_stove = 0;
    F_stove_shroud = 1 - F_stove_pot;
    
    %pot view factors
    F_pot_pot = 0;
    F_pot_stove = (A_stove / A_pot) * F_stove_pot;
    F_pot_shroud = 1 - F_pot_stove;
    
    %side view factors
    F_shroud_pot = (A_pot / A_shroud) * F_pot_shroud;
    F_shroud_stove = (A_stove / A_shroud) * F_stove_shroud;
    
    
    % Setup System for Finding q_pot
    sympref('FloatingPointOutput', true);
    
    sigma = 5.67*10^(-8);
    
    syms J_1 J_2 J_3;
    stove1 = J_1 == sigma*T_stove^4;
    pot2 = (sigma*T_pot^4 - J_2)*(epsilon_pot * A_pot) / (1-epsilon_pot) == A_pot * (F_pot_stove*(J_2-J_1) + F_pot_shroud*(J_2-J_3));
    shroud3 = (sigma*T_shroud^4 - J_3)*(epsilon_shroud * A_shroud) / (1-epsilon_shroud) == A_shroud * (F_shroud_stove*(J_3-J_1) + F_shroud_pot*(J_3-J_2));
    soln = solve([stove1 pot2 shroud3], [J_1 J_2 J_3]);

   if H > .099 && H < .101
       disp("Scenario 3 Radioosity Values");
       disp([soln.J_1 soln.J_2 soln.J_3]);
    end
    
    q_pot = (sigma*T_pot^4 - soln.J_2)*(epsilon_pot * A_pot) / (1-epsilon_pot);
    heatTransferArray3 = [heatTransferArray3, q_pot];
end

%% Plot Data

%plot view factor as a function of height
figure('Name', 'View Factor Stove to Pot v. Height')
hold on

plot(heightArray, viewFactorArray);

xlabel('Height [m]');
ylabel('View Factor');
title('View Factor Stove to Pot v. Height');

hold off

%plot heat transfer as function of height for different shrouds
figure('Name', 'Heat Transfer into Pot v. Height') 
hold on

%plot data, note the arrays are multiplied by -1 to represent the heat into
%the pot as a positive number
plot(heightArray, -1.*heatTransferArray1);
plot(heightArray, -1.*heatTransferArray2);
plot(heightArray, -1.*heatTransferArray3);
legend('No Shroud (virtual surface)', 'Aluminum Shroud (epsilon = .1)', 'Painted Shroud (epsilon = .9)');

xlabel('Height [m]');
ylabel('Heat Transfer [W]');
title('Heat Transfer into Pot v. Height');



hold off

