clc, clear, close all

Np = 10000; % number of particles
Nt = 500; %2100; % number of time steps
T = 1; % end time
nindx = [1,20,100,500]+1; % plot solution at t = (nindx-1)*tau
xint = 0.5; % interface location (must coincide with a lattice site)
tau = T/Nt; % time step

time_discretisation = 'exact exponential';
% time_discretisation = 'forward euler';

% Nr = number of stochastic model simulations/realisations
% stochastic_algorithm 1: tracks number of particles at each lattice site only
% stochastic_algorithm 2: tracks positions of individual particles

% animation = true; Nr = 1; stochastic_algorithm = '1';
animation = false; Nr = 5; stochastic_algorithm = '2';
colors = [128,0,0; 245,130,48; 0 0 128; 0 130 200]/255;

%% Problems
% L is length of domain and Nx is number of lattice sites
% model = '1'; L = 1; Nx = 101; xmin = 0.4; xmax = 0.6; D = 0.1; params.D = D;
% model = '2'; L = 1; Nx = 101; xmin = 0.4; xmax = 0.6; D1 = 0.1; D2 = 0.01; params.D1 = D1; params.D2 = D2;
% model = '3'; L = 1; Nx = 101; xmin = 0.2; xmax = 0.4; D1 = 0.1; D2 = 0.01; H = 0.5; params.D1 = D1; params.D2 = D2; params.H = H;
% model = '4'; L = 5; Nx = 501; xmin = 0.2; xmax = 0.4; D = 0.1; v = 1.0; params.D = D; params.v = v;
% model = '5'; L = 5; Nx = 501; xmin = 0.2; xmax = 0.4; D1 = 0.1; D2 = 0.01; v1 = 1.0; v2 = 1.0; params.D1 = D1; params.D2 = D2; params.v1 = v1; params.v2 = v2;
model = '6'; L = 5; Nx = 501; xmin = 0.2; xmax = 0.4; D1 = 0.1; D2 = 0.01; v1 = 1.0; v2 = 1.0; H = 0.5; params.D1 = D1; params.D2 = D2; params.v1 = v1; params.v2 = v2; params.H = H;

x = linspace(0,L,Nx)'; % node positions
f = @(x) 1.0*(x >= xmin & x <= xmax) + 0.0; % initial continuum particle density
delta = L/(Nx-1); % lattice spacing

if ismember(model,{'2','3','5','6'}) && isempty(find(abs(x-xint)/xint<1e-6,1,'first'))
    warning('Note delta must divide evenly into xint so xint/delta is an integer.')
end

% Impefect contact only
if isequal(model,'3') || isequal(model,'6')
    x = sort([x; xint]); % add in second lattice site at interface
    Ns = Nx + 1; % number of lattice sites
else
    Ns = Nx; % number of lattice sites
end

%% Discretisation in space using finite volume method
symbolic = false;
[A,V] = spatial_discretisation(x,Ns,delta,xint,model,params,symbolic);

%% Transition matrix
I = eye(Ns,Ns);
Vm = diag(V);
Vminv = diag(1./V);
C = (Vm*(A*Vminv))';
if isequal(time_discretisation,'exact exponential')
    P = expm(tau*C);
    minP = min(min(P));
elseif isequal(time_discretisation,'forward euler')
    P = I + tau*C;
end

%% Check constraints
if isequal(model,'1')
    max_tau = delta^2/(2*D);
elseif isequal(model,'2')
    max_tau = min(delta^2/(2*D1),delta^2/(2*D2));
elseif isequal(model,'3')
    max_tau = min(delta^2/(2*(D1+H*delta)),delta^2/(2*(D2+H*delta)));
elseif isequal(model,'4')
    if isequal(time_discretisation,'forward euler')
        max_tau = min(delta^2/(2*D+v*delta),delta^2/(2*D-v*delta));
    end
    max_delta = 2*D/abs(v);
elseif isequal(model,'5')
    if isequal(time_discretisation,'forward euler')
        max_tau = min([delta^2/(2*D1),delta^2/(2*D2),delta^2/(2*D1+v1*delta),delta^2/(2*D2-v2*delta),2*delta^2/(2*(D1+D2)+(v2-v1)*delta)]);
    end
    max_delta = min(2*D1/abs(v1),2*D2/abs(v2));
elseif isequal(model,'6')
    if isequal(time_discretisation,'forward euler')
        max_tau = min([delta^2/(2*D1),delta^2/(2*D2),delta^2/(2*D1+v1*delta),delta^2/(2*D2-v2*delta),delta^2/(2*D1-v1*delta+2*H*delta),delta^2/(2*D2+v2*delta+2*H*delta)]);
    end
    max_delta = min(2*D1/abs(v1),2*D2/abs(v2));
end

if ismember(model,{'1','2','3'}) && isequal(time_discretisation,'forward euler')
    if tau > max_tau
        warning(['To ensure non-negative probabilities for this model and choice of parameters, ',...
            'must have either tau <= %e (or Nt >= %i).'],max_tau,ceil(T/max_tau));
    end
elseif ismember(model,{'4','5','6'}) && isequal(time_discretisation,'forward euler')
    if tau > max_tau || delta > max_delta
        warning(['To ensure non-negative probabilities for this model and choice of parameters, ',...
            'must have both tau <= %e (or Nt >= %i) and delta <= %e (or Nx >= %i).'],max_tau,ceil(T/max_tau),...
            max_delta,ceil(L/max_delta+1));
    end
elseif ismember(model,{'4','5','6'}) && isequal(time_discretisation,'exact exponential')
    if delta > max_delta
        warning(['To ensure non-negative probabilities for this model and choice of parameters, ',...
            'must have delta <= %e (or Nx >= %i). \nDisclaimer: these bounds are approximate as ',...
            'the matrix exponential is computed numerically.'],max_delta,ceil(L/max_delta+1));
    end
end

%% Check transition matrix is a (right) stochastic matrix
fprintf('Each row of P sums to one: ');
if sum(abs(sum(P,2)-ones(Ns,1))<=1e-12)==Ns % check rows of P sum to one
    fprintf('YES\n');
else
    fprintf('NO\n');
end
fprintf('Each entry of P is non-negative: ');
if sum(sum(P>=0))==Ns^2 % check entries of P are all non-negative
    fprintf('YES\n');
else
    fprintf('NO\n');
end

%% Initial condition
Uc0 = f(x); % initial continuum particle density
Sp = (Uc0'*V)/Np; % scaling constant
N0 = round((Uc0.*V)/Sp); % initial number of particles at lattice sites
Np = sum(N0); % revised number of particles

N = zeros(Ns,Nt+1);
N(:,1) = N0;

stochastic_algorithm1 = isequal(stochastic_algorithm,'1');
stochastic_algorithm2 = isequal(stochastic_algorithm,'2');

if stochastic_algorithm2
    % Initial particle positions
    x0 = zeros(Np,1);
    for i = 1:Ns
        if i == 1
            x0(1:N0(1)) = i;
        else
            x0(sum(N0(1:i-1))+1:sum(N0(1:i-1))+N0(i)) = i;
        end
    end
    x0 = x0(randperm(Np));
    xp = zeros(Np,Nt+1);
    xp(:,1) = x0; % store lattice site number rather than position
    xp_store = zeros(Np,length(nindx),Nr);
end

Us_store = zeros(Ns,length(nindx),Nr);
Uc_store = zeros(Ns,length(nindx),Nr);

% cumulative probabilities
CP = zeros(size(P));
for i = 1:Ns
    CP(i,:) = cumsum(P(i,:));
end

expA = expm(tau*A); % matrix exponential for continuum particle density

%% Stochastic algorithm
for m = 1:Nr
    Uc = Uc0; % initial particle density (continuum)
    for n = 2:Nt+1

        if stochastic_algorithm1
            N(:,n) = N(:,n-1);
            for i = 1:Ns
                for k = 1:N(i,n-1)
                    r = rand;
                    j = find((CP(i,:)-r)>0,1,'first');
                    N(i,n) = N(i,n) - 1; N(j,n) = N(j,n) + 1;
                end
            end
        elseif stochastic_algorithm2
            for k = 1:Np
                r = rand;
                i = xp(k,n-1);
                j = find((CP(i,:)-r)>0,1,'first');
                xp(k,n) = j;
            end
            for i = 1:Ns
                N(i,n) = length(find(xp(:,n)==i));
            end
        end

        Us = N(:,n)*Sp./V; % stochastic particle density
        Uc = expA*Uc; % continuum particle density

        if animation
            % Plot stochastic and continuum models at each time step
            plot(x,Us,'Color',colors(1,:),'LineStyle','-','LineWidth',2)
            hold on
            plot(x,Uc,'Color',colors(3,:),'LineStyle','-','LineWidth',2);
            set(gca,'FontSize',16,'Xtick',[0,0.5,1],'YTick',[0,0.5,1])
            xlabel('Position ','FontSize',22)
            ylabel('Particle Density','FontSize',22)
            legend('Stochastic','Continuum')
            ylim([0,1.1])
            xlim([0,1])
            hold off
            drawnow
        else
            % Store solutions for specified times
            if ismember(n,nindx)
                [~,cnt] = ismember(n,nindx);
                Us_store(:,cnt,m) = Us;
                Uc_store(:,cnt,m) = Uc;
                if stochastic_algorithm2
                    xp_store(:,cnt,m) = xp(:,n);
                end
            end
        end
    end
    fprintf('Stochastic simulation %i of %i completed.\n',m,Nr);
end

%% Plot individual particle positions for one stochastic model simulation
if stochastic_algorithm2 && ~animation
    figure;
    fig1 = tiledlayout(4,1);
    for i = 1:length(nindx)
        nexttile;
        plot((xp_store(:,i,1)-1)*delta,1:Np,'.','Color',colors(i,:))
        set(gca,'FontSize',16,'Xtick',[0,0.5,1])
        hold on
        xlim([0,1])
        xlabel(fig1,'Position','FontSize',22)
        ylabel(fig1,'Particle Number','FontSize',22)
    end
    path_name = '/Users/carre/Dropbox/Documents/Research/Projects/1D Random Walk/Paper/Figures/';
    print(gcf,[path_name,'Problem_Stochastic',num2str(model)],'-depsc2')
end

%% Plot particle density from continuum model and stochastic model
if ~animation
    figure;
    leg = zeros(length(nindx),1); leglabels = cell(length(nindx),1);
    for i = 1:length(nindx)
        q1 = quantile(Us_store(:,i,:),0.025,3);
        q2 = quantile(Us_store(:,i,:),0.975,3);
        op_color = ones(1,3)-0.8*(ones(1,3)-colors(i,:));
        plot(x,q1,'Color',op_color,'LineStyle','-','LineWidth',1)
        hold on
        plot(x,q2,'Color',op_color,'LineStyle','-','LineWidth',1);
        patch([x,fliplr(x)]',[q1,fliplr(q2)]',1,'FaceColor',op_color,'EdgeColor',op_color);
        leg(i) = plot(x,Uc_store(:,i),'Color',colors(i,:),'LineStyle','-','LineWidth',2);
        set(gca,'FontSize',16,'Xtick',[0,0.5,1],'YTick',[0,0.5,1])
        xlabel('Position ','FontSize',22)
        ylabel('Particle Density','FontSize',22)
        ylim([0,1.1])
        xlim([0,1])
        box on
        drawnow
        leglabels{i} = [' {\itt} = ',num2str((nindx(i)-1)*tau,'%g')];
    end
    legend(leg,leglabels)
    path_name = '/Users/carre/Dropbox/Documents/Research/Projects/1D Random Walk/Paper/Figures/';
    print(gcf,[path_name,'Problem',num2str(model)],'-depsc2')
end