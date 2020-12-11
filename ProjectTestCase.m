% This is a test case that I used to troubleshoot this code. It also
% shows what a large deformation would look like (hard to tell if the
% deformed shape plot is working with the small deformations that occur 
% in the 5 provided cases).

clear all; close all; clc; % Clear everything, close everything.
% TRUSS PARAMETERS ENTERED HERE
nodes = [0 0;4 0;8 0;12 0;16 0;20 0;24 0;4 4;20 4;8 8;12 8;16 8]; % Enter the coordinates of our joints (meters) [x y].
elements = [1 2;1 8;2 3;2 8;3 4;3 10;4 5;4 12;4 11;5 6;5 9;5 12;6 7;6 9;8 3;8 10;9 7;10 4;10 11;11 12;12 9]; % Enter the joints that our members connect [startjoint endjoint].
E = 45000000000; % Enter the Modulus of Elasticity here (Pa). For magnesium alloy, this is 45 GPa.
A = 0.00755; % Enter our area here (m^2).
sigmay = 130000000; % Enter our yield strength here (for magnesium alloy, this is 130 MPa for both compressive and tensile stresses).
unconstrained = [3 4 5 6 7 8 9 10 11 12 15 16 17 18 19 20 21 22 23 24]; % Enter boundary conditions: the unconstrained joint directions (3 4 means joint 2 is free in x and y).
constrained = [1 2 13 14]; % Enter boundary conditions: the constrained joint directions (1 2 means joint 1 cannot move in x and y).
% CTest
load1 = [0;0;0;0;0;0;0;10000000;0;10000000;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
b1=-9.81*load1;

figure1=figure('Position', [100, 100, 800, 600]); % Plot options including vertical/horizontal size.
[maxstress(1,1)] = TrussDirectStiffness(nodes,elements,b1,E,A,unconstrained,constrained,1,sigmay,1);