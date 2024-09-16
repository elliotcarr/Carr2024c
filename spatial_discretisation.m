function [A,V] = spatial_discretisation(x,Ns,delta,xint,model,params,symbolic)
% Inputs
%    D: Diffusivity (constant scalar)
%    Ns: Number of nodes (positive integer) (Assume Ns >= 3, so you have at least one interior node)
%    t: Discrete times (vector)
%    u0: Initial numerical solution (vector)
% Outputs
%    u: Solution (matrix of size Ns by length(t)), where u(i,j) is the solution u(x,t) evaluated
%       at x = (i-1)h and t = t(j) with h = 1/(Ns-1)

if symbolic
    A = sym(zeros(Ns,Ns));
    V = sym(zeros(Ns,1));
    Nx1 = sym(find(x==xint,1,'first')); % first interface lattice site
else
    A = zeros(Ns,Ns);
    V = zeros(Ns,1);
    Nx1 = find(abs(x-xint)/xint<1e-6,1,'first'); % first interface lattice site
end

% Control volume lengths
V(1) = delta/2; V(2:Ns-1) = delta; V(Ns) = delta/2;

if isequal(model,'3') || isequal(model,'6') 
    V(Nx1) = delta/2; V(Nx1+1) = delta/2;
end

%% Model 1: Homogeneous diffusion
if isequal(model,'1')
    D = params.D;
    A(1,1) = -2*D/delta^2;
    A(1,2) = 2*D/delta^2;
    for i = 2:Ns-1
        A(i,i-1) = D/delta^2;
        A(i,i) = -2*D/delta^2;
        A(i,i+1) = D/delta^2;
    end
    A(Ns,Ns-1) = 2*D/delta^2;
    A(Ns,Ns) = -2*D/delta^2;
end

%% Model 2: Heterogeneous diffusion with fully-permeable interface
if isequal(model,'2')
    D1 = params.D1;
    D2 = params.D2;
    A(1,1) = -2*D1/delta^2; 
    A(1,2) = 2*D1/delta^2;
    for i = 2:Nx1-1
        A(i,i-1) = D1/delta^2;
        A(i,i) = -2*D1/delta^2;
        A(i,i+1) = D1/delta^2;
    end
    for i = Nx1
        A(i,i-1) = D1/delta^2;
        A(i,i) = -(D1+D2)/delta^2;
        A(i,i+1) = D2/delta^2;
    end
    for i = Nx1+1:Ns-1
        A(i,i-1) = D2/delta^2;
        A(i,i) = -2*D2/delta^2;
        A(i,i+1) = D2/delta^2;
    end
    A(Ns,Ns-1) = 2*D2/delta^2; 
    A(Ns,Ns) = -2*D2/delta^2;
end

%% Model 3: Heterogeneous diffusion with semi-permeable interface
if isequal(model,'3')
    D1 = params.D1;
    D2 = params.D2;
    H = params.H;    
    A(1,1) = -2*D1/delta^2;
    A(1,2) = 2*D1/delta^2;
    for i = 2:Nx1-1
        A(i,i-1) = D1/delta^2;
        A(i,i) = -2*D1/delta^2;
        A(i,i+1) = D1/delta^2;
    end
    for i = Nx1
        A(i,i-1) = 2*D1/delta^2;
        A(i,i) = -2*(D1 + H*delta)/delta^2;
        A(i,i+1) = 2*H/delta;
    end
    for i = Nx1+1
        A(i,i-1) = 2*H/delta;
        A(i,i) = -2*(D2 + H*delta)/delta^2;
        A(i,i+1) = 2*D2/delta^2;
    end
    for i = Nx1+2:Ns-1
        A(i,i-1) = D2/delta^2;
        A(i,i) = -2*D2/delta^2;
        A(i,i+1) = D2/delta^2;
    end
    A(Ns,Ns-1) = 2*D2/delta^2; 
    A(Ns,Ns) = -2*D2/delta^2;
end

%% Model 4: Homogeneous advection-diffusion
if isequal(model,'4')
    D = params.D;
    v = params.v;
    A(1,1) = -(2*D + v*delta)/delta^2; 
    A(1,2) = (2*D - v*delta)/delta^2;
    for i = 2:Ns-1
        A(i,i-1) = (2*D + v*delta)/(2*delta^2);
        A(i,i) = -2*D/delta^2;
        A(i,i+1) = (2*D - v*delta)/(2*delta^2);
    end
    A(Ns,Ns-1) = (2*D + v*delta)/delta^2;
    A(Ns,Ns) = -(2*D - v*delta)/delta^2;
end

%% Model 5: Heterogeneous advection-diffusion with fully-permeable interface
if isequal(model,'5')
    D1 = params.D1;
    D2 = params.D2;
    v1 = params.v1;
    v2 = params.v2;
    A(1,1) = -(2*D1 + v1*delta)/delta^2; 
    A(1,2) = (2*D1 - v1*delta)/delta^2;    
    for i = 2:Nx1-1
        A(i,i-1) = (2*D1 + v1*delta)/(2*delta^2);
        A(i,i) = -2*D1/delta^2;
        A(i,i+1) = (2*D1 - v1*delta)/(2*delta^2);        
    end
    for i = Nx1
        A(i,i-1) = (2*D1 + v1*delta)/(2*delta^2);
        A(i,i) = -(2*(D1 + D2) + (v2 - v1)*delta)/(2*delta^2);
        A(i,i+1) = (2*D2 - v2*delta)/(2*delta^2);        
    end
    for i = Nx1+1:Ns-1
        A(i,i-1) = (2*D2 + v2*delta)/(2*delta^2);
        A(i,i) = -2*D2/delta^2;
        A(i,i+1) = (2*D2 - v2*delta)/(2*delta^2);        
    end
    A(Ns,Ns-1) = (2*D2 + v2*delta)/delta^2; 
    A(Ns,Ns) = -(2*D2 - v2*delta)/delta^2;    
end

%% Model 6: Heterogeneous advection-diffusion with semi-permeable interface
if isequal(model,'6')
    D1 = params.D1;
    D2 = params.D2;
    v1 = params.v1;
    v2 = params.v2;
    H = params.H;
    A(1,1) = -(2*D1 + v1*delta)/delta^2; 
    A(1,2) = (2*D1 - v1*delta)/delta^2;    
    for i = 2:Nx1-1
        A(i,i-1) = (2*D1 + v1*delta)/(2*delta^2);
        A(i,i) = -2*D1/delta^2;
        A(i,i+1) = (2*D1 - v1*delta)/(2*delta^2);        
    end
    for i = Nx1       
        A(i,i-1) = (2*D1 + v1*delta)/delta^2;
        A(i,i) = -(2*D1 - v1*delta + 2*H*delta)/delta^2;
        A(i,i+1) = 2*H/delta;
    end
    for i = Nx1+1
        A(i,i-1) = 2*H/delta;
        A(i,i) = -(2*D2 + v2*delta + 2*H*delta)/delta^2;
        A(i,i+1) = (2*D2 - v2*delta)/delta^2;        
    end
    for i = Nx1+2:Ns-1
        A(i,i-1) = (2*D2 + v2*delta)/(2*delta^2);
        A(i,i) = -2*D2/delta^2;
        A(i,i+1) = (2*D2 - v2*delta)/(2*delta^2);    
    end
    A(Ns,Ns-1) = (2*D2 + v2*delta)/delta^2; 
    A(Ns,Ns) = -(2*D2 - v2*delta)/delta^2;    
end