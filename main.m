% Heat Transfer - Spring 2023 Computational Assignment 3
% Michael Allen, Cullen Hirstius, Hayden Payne
% 4/21/23

%% Scenario 1
%loop through all values of H once corrected

%% Problem Setup

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
H = .1; %height from flame to pot [m]
A_surr = pi*diameter_pot * H;


%% Calculate View Factors
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


%% Setup System for Finding q_pot
sympref('FloatingPointOutput', true);

sigma = 5.67*10^(-8);

syms J_1 J_2 J_3;
stove1 = J_1 == sigma*T_stove^4;
pot2 = (sigma*T_pot^4 - J_2)*(epsilon_pot * A_pot) / (1-A_pot) == A_pot * (F_pot_stove*(J_2-J_1) + F_pot_surr*(J_2-J_3));
surroundings3 = J_3 == (sigma)*T_surr^4;
soln = solve([stove1 pot2 surroundings3], [J_1 J_2 J_3]);

q_pot = (sigma*T_pot^4 - soln.J_2)*(epsilon_pot * A_pot) / (1-A_pot);
disp(q_pot);

%% Scenario 2