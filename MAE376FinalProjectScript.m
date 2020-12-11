% An analysis of nodal reaction forces, nodal deformations, element
% internal forces, and element internal stresses for a truss of a certain
% configuration is required. This analysis must use the direct stiffness
% method. A plot of maximum stress vs loaded mass is also desired, as is a
% deformed shape plot for Case 5.

% For this particular truss (UBID 50084114), height and width are both 
% 4 meters. It is composed of magnesium alloy, and the masses loaded
% are 150 kg. Material properties are assumed to be those of Magnesium
% AM60B Cast Alloy (the properties of this alloy can be found at the
% following link: http://www.azom.com/article.aspx?ArticleID=9237#4). Case
% 1 has a mass at joint 2. Case 2 has a mass at joint 2 and 3. Case 3 has a
% mass at joint 2, 3, and 4. Case 4 has a mass at joint 2, 3, 4, and 5.
% Case 5 has a mass at joint 2, 3, 4, 5, and 6.

% The effects of member weight and member buckling are ignored in this
% analysis, and all joints are assumed to be pin joints.

% The truss was analyzed using a script and function. The script stores
% user inputs (load vector, E, A, nodes, etc) for the parameters of the
% truss, and passes them along to the function TrussDirectStiffness. It is
% also responsible for plotting maximum stress vs mass loaded, and
% creating a fit function along with calculating the goodness of fit for
% that function.

% The TrussDirectStiffness function takes user inputs for the b vector,
% along with the node coordinates, members, E, A, and other information,
% and uses this to create localized stiffness matrices and concatenate
% them. Finally, it solves the resulting system of equations (using
% boundary conditions to eliminate rows and columns of the K matrix and b
% and u vectors that would otherwise render the system unsolvable) and
% returns a variety of outputs.

% For each case, the reactions, nodal deformations, internal forces, and 
% internal stresses are outputted to the user. The highest magnitude axial 
% stresses and nodal displacements are also shown. A yielding check is
% performed (if the user so desires) as a preliminary check for member
% failure. The MATLAB function TrussDirectStiffness also outputs a shape 
% plot if prompted to do so, which shows a representation of the truss in its
% deformed state.

% For the analysis results of each case, along with a test case showing the
% deformed shape plot in greater detail, see the code output. Design
% details are included in the extensive comments in the code itself.

clear all; close all; clc; % Clear everything, close everything.
% TRUSS PARAMETERS ENTERED HERE
nodes = [0 0;4 0;8 0;12 0;16 0;20 0;24 0;4 4;20 4;8 8;12 8;16 8]; % Enter the coordinates of our joints (meters) [x y].
elements = [1 2;1 8;2 3;2 8;3 4;3 10;4 5;4 12;4 11;5 6;5 9;5 12;6 7;6 9;8 3;8 10;9 7;10 4;10 11;11 12;12 9]; % Enter the joints that our members connect [startjoint endjoint].
E = 45000000000; % Enter the Modulus of Elasticity here (Pa). For magnesium alloy, this is 45 GPa.
A = 0.00755; % Enter our area here (m^2).
sigmay = 130000000; % Enter our yield strength here (for magnesium alloy, this is 130 MPa for both compressive and tensile stresses).
unconstrained = [3 4 5 6 7 8 9 10 11 12 15 16 17 18 19 20 21 22 23 24]; % Enter boundary conditions: the unconstrained joint directions (3 4 means joint 2 is free in x and y).
constrained = [1 2 13 14]; % Enter boundary conditions: the constrained joint directions (1 2 means joint 1 cannot move in x and y).
% END OF TRUSS PARAMETERS
% C1
load1 = [0;0;0;150;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; % Enter the load vector (kg).
b1=-9.81*load1; % The load vector is converted into N, and comprises the initial b vector fed into the function.
[maxstress(1,1)] = TrussDirectStiffness(nodes,elements,b1,E,A,unconstrained,constrained,1,sigmay); % Call TrussDirectStiffness function with defined parameters, and obtain the max stress (which is stored in a column vector).
% C2
load2 = [0;0;0;150;0;150;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
b2=-9.81*load2;
[maxstress(2,1)] = TrussDirectStiffness(nodes,elements,b2,E,A,unconstrained,constrained,1,sigmay);
% C3
load3 = [0;0;0;150;0;150;0;150;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
b3=-9.81*load3;
[maxstress(3,1)] = TrussDirectStiffness(nodes,elements,b3,E,A,unconstrained,constrained,1,sigmay);
% C4
load4 = [0;0;0;150;0;150;0;150;0;150;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
b4=-9.81*load4;
[maxstress(4,1)] = TrussDirectStiffness(nodes,elements,b4,E,A,unconstrained,constrained,1,sigmay);
% C5
load5 = [0;0;0;150;0;150;0;150;0;150;0;150;0;0;0;0;0;0;0;0;0;0;0;0];
b5=-9.81*load5;

figure1=figure('Position', [100, 100, 800, 600]); % Plot options including vertical/horizontal size.
[maxstress(5,1)] = TrussDirectStiffness(nodes,elements,b5,E,A,unconstrained,constrained,1,sigmay,1); % Note that the deformed shape plot has been enabled.

mass(1,1) = sum(load1); % Calculate the total weight loaded on the truss for each loading scenario.
mass(2,1) = sum(load2);
mass(3,1) = sum(load3);
mass(4,1) = sum(load4);
mass(5,1) = sum(load5);

figure2=figure('Position', [100, 100, 800, 600]); % Plot options including vertical/horizontal size.
scatter(mass,abs(maxstress),'filled') % A scatter plot is most appropriate for these discrete points, this one plots the max stress vs total mass loaded.
xlabel('Total Mass (kg)','fontsize',16);
ylabel('Maximum Stress Magnitude (Pa)','fontsize',16);
title('Maximum Stress Magnitude vs Total Mass','fontsize',16); % Title, x label, y label, x interval settings, etc.
set(gca,'XTick',0:150:mass(end,1));
hold on

p = polyfit(mass(:,1),abs(maxstress(:,1)),3); % Fit the max stress vs mass scatterplot with a polynomial, order 3.
yfit = @(x) p(1)*x.^3 + p(2)*x.^2 + p(3)*x +p(4); % Using the four coefficients outputted by the built in polyfit function, we define fit function.
fplot(yfit,[mass(1,1) mass(end,1)]); % Plot our fit function.

resid = abs(maxstress(:,1)) - feval(yfit,mass(:,1)); % This series of calculations evaluates goodness of fit.
SSresid = sum(resid.^2);
SStotal = (length(maxstress(:,1))-1)*var(maxstress(:,1));
r2 = 1-SSresid/SStotal; % r2, the coefficient of determination, is calculated.

disp(['The coefficient of determination r^2 is ' num2str(r2) ' for this fitting function.']); % Text output for user.

