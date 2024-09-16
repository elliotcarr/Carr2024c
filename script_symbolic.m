clc, clear, close all

syms delta real positive
Nx = 7; % number of lattice sites
x = (0:Nx-1)'*delta; % lattice site positions
xint = 3*delta; % interface location (must coincide with a lattice site)

% model = '1'; 
% model = '2'; 
% model = '3';
% model = '4'; 
% model = '5'; 
model = '6'; 

if isequal(model,'1')
    syms D real positive
    syms tau real positive
    params.D = D;
elseif isequal(model,'2')
    syms D1 real positive
    syms D2 real positive
    syms tau real positive
    params.D1 = D1; params.D2 = D2;
elseif isequal(model,'3')
    syms D1 real positive
    syms D2 real positive
    syms H real positive
    syms tau real positive
    params.D1 = D1; params.D2 = D2; params.H = H;
elseif isequal(model,'4')
    syms D real positive
    syms v real
    syms tau real positive
    params.D = D; params.v = v;
elseif isequal(model,'5')
    syms D1 real positive
    syms D2 real positive
    syms v1 real positive
    syms v2 real positive
    syms tau real positive
    params.D1 = D1; params.D2 = D2; params.v1 = v1; params.v2 = v2;
elseif isequal(model,'6')
    syms D1 real positive
    syms D2 real positive
    syms v1 real positive
    syms v2 real positive
    syms H real positive
    syms tau real positive
    params.D1 = D1; params.D2 = D2; params.v1 = v1; params.v2 = v2; params.H = H;   
end

% Impefect contact only
if isequal(model,'3') || isequal(model,'6')
    x = sort([x; xint]); % add in second lattice site at interface
    Ns = Nx + 1; % number of lattice sites
else 
    Ns = Nx; % number of lattice sites
end

symbolic = true;
[A,V] = spatial_discretisation(x,Ns,delta,xint,model,params,symbolic);
display(A)

%% Discretisation matrix
I = eye(Ns,Ns);
Vm = diag(V);
Vminv = diag(1./V);
C = (Vm*(A*Vminv))';
display(C)

%% Transition matrix (forward Euler discretisation)
P = I + tau*C;
P = expand(P);
display(P)
fprintf('sum(P,2)\n')
simplify(sum(P,2))
fprintf('sum(C,2)\n')
simplify(sum(C,2))
fprintf('isequal(A,C)\n')
isequal(A,C)