function [maxstress maxdisplacement reaction displacement force stress b u f st K] = TrussDirectStiffness(nodes,elements,b,E,A,unconstrained,constrained,yieldcheck,sigmay,shape) % Define our inputs, outputs.

% Function outputs max stress, max displacement, reaction matrix,
% displacement matrix, internal force matrix, axial stress matrix, as 
% well as raw form force vector, nodal displacement vector, internal 
% force vector, and axial stress vector. The global stiffness matrix K is
% also outputted for troubleshooting purposes. The function also performs a
% preliminary yield check for all members using inputted yield strength
% value (if one is inputted), and outputs a deformed shape plot if prompted
% to do so.

% Function requires nodes matrix, elements matrix, b vector, E, A,
% unconstrained and constrained matrix (boundary conditions for the direct
% stiffness method.

if nargin<7,error('The following input arguments are required: nodes,elements,b,E,A,unconstrained,constrained'),end % Check for sufficient input arguments.
if nargin==8 && yieldcheck ~= 0,error('To enable yield check, enter the yield strength of the material sigmay. Otherwise, yieldcheck = 0 or simply leaving the yieldcheck input blank will disable yield check.'),end % Enabling yield check without inputting sigma y causes an error.
if nargin<8,yieldcheck = 0; end % If there is no yieldcheck input provided, that loop is disabled.
if nargin<10,shape = 0; end % If the user does not request a deformed shape plot, disable loop.
if yieldcheck ~= 0 && yieldcheck ~= 1, error('Input for yieldcheck must be 0 (off) or 1 (on)'),end
if sigmay <= 0, error('Yield strength must be positive.'),end
if shape ~= 0 && shape ~= 1, error('Input for shape must be 0 (off) or 1 (on)'),end
if size(nodes,2) ~= 2, error('Node coordinate matrix must be of size n by 2.'),end % Next lines check for dimensional mismatches.
if size(elements,2) ~= 2, error('Element matrix must be of size n by 2.'),end
if size(b,2)~= 1, error('b must be a column vector.'),end
if size(b,1)~= 2*size(nodes,1), error('b must be a column vector of size 2n'),end
if size(unconstrained,1)~= 1, error('unconstrained must be a row vector.'),end
if size(constrained,1)~= 1, error('constrained must be a row vector.'),end 
if size(unconstrained,2)+size(constrained,2) ~= 2*size(nodes,1), error('Not enough boundary conditions specified.'),end
if size(E,1)~= 1, error('E must be a scalar value.'),end
if size(E,2)~= 1, error('E must be a scalar value.'),end
if size(A,1)~= 1, error('A must be a scalar value.'),end
if size(A,2)~= 1, error('A must be a scalar value.'),end

Nelements = size(elements,1); % The number of "elements" for the direct stiffness method (an element is a truss member linking two joints).
Nnodes = size(nodes,1); % The number of "nodes" for the direct stiffness method (a node is a truss joint).

K = zeros(2*Nnodes); % Initializing K matrix, the assembled global stiffness matrix (must be square matrix of size = 2*number of nodes, because each node has vertical and horizontal degree of freedom). Values are unknown, therefore initialize as zeros.
u = zeros(2*Nnodes,1); % Initializing u displacement matrix (column vector of size = 2*number of nodes, 1). Values are unknown, therefore initialize as zeros.

for i=1:Nelements % Iterate from 1 to the number of elements.
    elementnodes = elements(i,1:2); % This selects the start and end nodes (start and end joints) of each member.
    nodecoordinates = nodes(elementnodes,:); % This selects the corresponding node coordinates from the joints matrix.
    
    x1 = nodecoordinates(1,1); % Retrieves the x1, x2, y1, y2 values from the supplied nodecoordinates matrix.
    x2 = nodecoordinates(2,1);
    y1 = nodecoordinates(1,2);
    y2 = nodecoordinates(2,2);

    L = sqrt((x2-x1)^2+(y2-y1)^2);  % Length is calculated in terms of node coordinates.
    cos = (x2-x1)/L; % Calculation of cos and sin for usage in Klocal matrix.
    sin = (y2-y1)/L; 

    Klocal = (E*A/L)*[cos^2 cos*sin -cos^2 -cos*sin;cos*sin sin^2 -cos*sin -sin^2;-cos^2 -cos*sin cos^2 cos*sin;-cos*sin -sin^2 cos*sin sin^2]; % The Klocal matrix is the stiffness matrix for one bar type element in global coordinates.
    
    idBeg = 2*(elementnodes(1)-1)+1:2*(elementnodes(1)-1)+2; % The displacement corresponding to the beginning of the element (what position it occupies in the u column vector).
    idEnd = 2*(elementnodes(2)-1)+1:2*(elementnodes(2)-1)+2; % The displacement corresponding to the end of the element.
    id = [idBeg idEnd]; % This identifies which values of the displacement vector u this local stiffness matrix corresponds to, and will be used in assembling the matrix below.
    K(id,id) = K(id,id) + Klocal; % The elementdisplacement vector is used to assemble the local stiffness matrices into one global stiffness matrix.
end

u(unconstrained) = K(unconstrained,unconstrained)\(b(unconstrained)-K(unconstrained,constrained)*u(constrained)); % The unconstrained node elements of the u vector are found (these are our unknown deformations).
b(constrained) = K(constrained,:)*u; % The elements of the force vector corresponding to constrained nodes are found (these are our unknown reaction forces).

for i=1:Nelements
    elementnodes = elements(i,1:2); % This selects the start and end nodes (start and end joints) of each member.
    nodecoordinates = nodes(elementnodes,:); % This selects the corresponding node coordinates from the joints matrix.
    
    x1 = nodecoordinates(1,1); % Retrieves the x1, x2, y1, y2 values from the supplied nodecoordinates matrix
    x2 = nodecoordinates(2,1);
    y1 = nodecoordinates(1,2);
    y2 = nodecoordinates(2,2);

    L = sqrt((x2-x1)^2+(y2-y1)^2);  % Length is calculated in terms of node coordinates.
    cos = (x2-x1)/L; % Calculation of cos and sin for usage in Klocal matrix.
    sin = (y2-y1)/L; 

    idBeg = 2*(elementnodes(1)-1)+1:2*(elementnodes(1)-1)+2; % The displacement corresponding to the beginning of the element (what position it occupies in the u column vector).
    idEnd = 2*(elementnodes(2)-1)+1:2*(elementnodes(2)-1)+2; % The displacement corresponding to the end of the element.
    
    f(i,1)=(E*A/L)*(cos*(u(idEnd(1),1)-u(idBeg(1),1))+sin*(u(idEnd(2),1)-u(idBeg(2),1))); % The internal force in each member is calculated based on the displacements of the beginning and end nodes, and arranged in a column vector.
    st(i,1) = (1/A)*f(i,1); % The stress is easily calculated, simply multiply the force vector by the scalar 1/A (equivalent to elementwise division by area, sigma = F/A).
end

[~,maxstresselement] = max(abs(st)); 
maxstress = st(sub2ind(size(st),maxstresselement,1:size(st,2))); % The greatest magnitude value in the matrix is returned, with +/- symbol

[~,maxdisplacementnode] = max(abs(u)); 
maxdisplacement = u(sub2ind(size(u),maxdisplacementnode,1:size(u,2))); % This process is repeated for maximum displacement.

reaction = zeros(Nnodes,3); % The u and b column vectors are broken down into matrix of size = (Nnodes,2) for readability, so that [x1 y1;x2 y2;...] corresponds to the x and y reactions/displacements of node 1, node 2, etc.
reaction(:,1) = 1:Nnodes;
reaction(:,2) = b(1:2:end);
reaction(:,3) = b(2:2:end);
displacement = zeros(Nnodes,3);
displacement(:,1) = 1:Nnodes;
displacement(:,2) = u(1:2:end);
displacement(:,3) = u(2:2:end);
force = zeros(Nelements,2); % The internal force and stress vectors are also placed into a more readable, indexed format.
force(:,1) = 1:Nelements;
force(:,2) = f(:,1);
stress = zeros(Nelements,2);
stress(:,1) = 1:Nelements;
stress(:,2) = st(:,1);

format short g % The reformatted reaction, displacement, force, and stress matrices are presented to the user.
disp(['Nodal Reactions (N) ----- [Node     X     Y]']); reaction
disp(['Nodal Displacements (m) ----- [Node     X     Y]']); displacement
disp(['Element Internal Forces (N) ----- [Element     Internal Force]']); force
disp(['Element Axial Stresses (Pa) ----- [Element     Axial Stress]']); stress
disp(['The largest magnitude axial stress is ' num2str(maxstress) ' (Pa).']); % Maximum stress (greatest magnitude) is displayed. Note that this may be compressive or tensile, and includes +/- symbol.
disp(['The largest magnitude displacement is ' num2str(maxdisplacement) ' (m).']); % Maximum displacement (greatest magnitude) is displayed.

if yieldcheck == 1 % If the user has decided to enable yieldchecking, this loop proceeds.
if abs(maxstress) > sigmay % This if loop checks maximum internal stress against the yield strength inputted into our function. 
    disp(['Yielding due to axial stress occurs in at least one member.']); % If yield strength is exceeded, the user is warned.
else
    disp(['Yielding due to axial stress does not occur.']); % If yield strength is not exceeded, no warning.
end
elseif yieldcheck == 0 % If the user has disabled yieldchecking or left yieldcheck input blank, the loop is disabled.
    return
end

if shape == 1 % If the "shape" input is 1, this triggers the loop that draws the deformed truss shape.
    deformedshape(:,1:2) = nodes(:,1:2)+ displacement(:,2:3); % The new positions of the nodes are calculated by adding the initial node positions and displacements.

    scatter(deformedshape(:,1),deformedshape(:,2),'filled') % A scatter plot is made of the nodes in their post deformation positions.
    xlabel('x (meters)','fontsize',16);
    ylabel('y (meters)','fontsize',16);
    title('Deformed Truss Shape','fontsize',16); % Plot settings and labels.
    hold on
    
for i=1:Nelements % Iterate from 1 to number of elements.
    elementnodes = elements(i,1:2); % Obtain the start and end nodes of each element (similar to the step that creates the global stiffness matrix).
    nodecoordinates = deformedshape(elementnodes,:); % This selects the corresponding node coordinates from the joints matrix.
    
    x1 = nodecoordinates(1,1); % Retrieves the x1, x2, y1, y2 values from the supplied nodecoordinates matrix
    x2 = nodecoordinates(2,1);
    y1 = nodecoordinates(1,2);
    y2 = nodecoordinates(2,2);
    
    m = (y2-y1)/(x2-x1); % This calculates the slope of the line formed by these nodes.
    y = @(x) m*x-m*x1+y1; % The equation of the line passing through both nodes is defined.
    if x2 > x1 
        fplot(y,[x1 x2]); % If x2>x1, the function is plotted on this range. Alternatively, the x2 and x1 positions are switched.
    elseif x2 < x1
        fplot(y,[x2 x1]);
    elseif x2 == x1
        if y2 > y1
        line([x1 x2],[y1 y2]); % If x1 and x2 are equal, a vertical line must be drawn. Note the y1 and y2 positions switch much like the x1 and x2 positions above.
        else
        line([x1 x2],[y2 y1]);
        end
    end
    set(gca,'XTick',min(deformedshape(:,1)):max(deformedshape(:,1))); % X and Y intervals are set.
    set(gca,'YTick',min(deformedshape(:,2)):max(deformedshape(:,2)));
    hold on
end
elseif shape == 0 % If the shape input is 0, the loop above is not triggered and a graph is not produced.
    return
end
end
