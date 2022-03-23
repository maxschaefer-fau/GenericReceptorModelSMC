% General Receptor model for GC 22 
%
% M. Schäfer, S. Lotter, M. Barros, 2022 
% Further details can be found in: 
% 
%
% General 3-state receptor model (Open, Closed, Desensitized) with
% possible saturation from O->C and O->D. The model can be specialized to
% any 3-state receptor by inserting the desired transition rates. 


%% Paramter Definition (with physical Dimensions) 

T_ = 1e-7;       % simulatiton time step in s 
t_end_ = 6e-3;   % simulation duration in s 
t_ = 0:T_:t_end_;  % time vector 

Na = 6.022e23; % Avogadro constant in mol^-1

a_ = 2*1e-2*1e-6; % channel width in x-direction in m
yz_ = 0.15;  % channel width in y- and z-direction in mum (1e-6 m) 

D_ = 3.3e-4*1e-12/1e-6; % Diffusion coefficient in m^2/s
keCe_ = 1e-3/1e-6; % Enzymatic degradation rate in 1/s

%% Receptor Definition (with physical Dimensions)
recepDensity_ = 1e3/1e-12;  % Receptor density in 1/m^2 TODO: Add reference!!
C = 203;                   % number of receptors 
r_rec_ = 2.3e-3;            % receptor radius in mum
rho_ = (C*pi*r_rec_^2)/yz_^2; % dimensionless receptor coverage


% (This is the AMPA setting from [...]) 
% C <---> O 
% kco_ = 0.995*1.02e-4*rho_;   % in m/s
% koc_ = 8.5e-3/1e-6;          % in 1/s
% %O <---> D 
% kod_ = 8.5e-3/1e-6;          % in 1/s
% kdo_ = 0;                    % in 1/s
% %C <---> D 
% kcd_ = 0;                    % in m/s
% kdc_ = 8.5e-3/1e-6;            % in 1/s

% NMDA Setting 
% C <---> O 
kco_ = 0;   % in m/s
koc_ = 8.5e-3/1e-6;          % in 1/s
% O <---> D 
kod_ = 0;          % in 1/s
kdo_ = 8.5e-3/1e-6;                    % in 1/s
% C <---> D 
kcd_ = 0.995*1.02e-4*rho_;                    % in m/s
kdc_ = 8.5e-3/1e-6;            % in 1/s

%% Normalization 
% Normalize all parameters to be dimensionless by Diffusion coefficient D
% and channel with in x-direction a

T = T_*D_/a_^2; 
t_end = t_end_*D_/a_^2; 
t = t_*D_/a_^2;

a = a_/a_; 

keCe = keCe_*a_^2/D_; 
D = D_/D_;

kco = kco_*a_/D_; 
koc = 0.5*koc_*a_^2/D_; 
kod = kod_*a_^2/D_; 
kdo = 0.1*kdo_*a_^2/D_; 
kcd = kcd_*a_/D_; 
kdc = 0.5*kdc_*a_^2/D_; 


%% Set up transfer function model for the inner synaptic model 

% see below eq. (26) in [...]
Mu = 20;           % number of eigenvalues 
mu = 0:Mu-1;        % vector to count eigenvalues 
gmu = mu*pi/a;      % wave numbers gamma 
smu = -D*gmu.^2;    % eigenvalues of the inner synaptic model 

% Eigenfunctions, see eq. (36) in [...]
K1 = @(x) cos(gmu*x);           % first entry of primal eigenfunction
K2 = @(x) D*gmu.*sin(gmu*x);    % second entry of primal eigenfunction

Ka1 = @(x) -D*gmu.*sin(gmu*x);  % first entry of adjoint eigenfunction
Ka2 = @(x) cos(gmu*x);          % second entry of adjoint eigenfunction

% Scaling factor, see eq. (27) 
nmu = ones(1,length(mu))*a/2;   % scaling factor due to bi-orthogonality of eigenfunctions 
nmu(1) = a;                     % scaling factor for mu = 0 has to be treated separately 

%% Set up matrices and vectors for TFM state space description 

% evaluate eigenfunctions at x = a (position where receptors are attached) 
ca1_a = Ka1(a); 
ca2_a = Ka2(a); 
c1_a = 1./nmu.*K1(a); 
c2_a = 1./nmu.*K2(a);

A = diag(exp(smu*T));   % matrix of eigenvalues in discrete-tie domain 
KeCe = exp(-keCe*T);    % enzymatic degradation --> exponential decay 


ybar = zeros(Mu,length(t));     % time evoluation of system states 
y = zeros(2,length(t));         % output vector y(t) = [c(t) i(t)]
ca = zeros(1,length(t));        % concentration at x = a, i.e., c(a,t)

% vectors for receptor fluxes and helping vectors 
io = zeros(1,length(t));        % number of receptors in open state 
id = zeros(1,length(t));        % number of receptors in desensitized state
ib = zeros(1,length(t));        % complete boundary flux 
phib = zeros(Mu,length(t));     % general boundary value 

iobar = zeros(1,length(t));     
idbar = zeros(1,length(t));

csat = zeros(1,length(t)); 
cs = zeros(1,length(t));

%% Excitation and inital values 
N0 = 1000;  % number of release NTs 

xe = 0;     % excitation position 
fx = N0*Ka2(xe); 

ft = zeros(1,length(t));    % temporal excitation function 
ft(1) = 1;                  % Define impulse at t = 0
ft(10001) = 1;                  % Define impulse at t = 1ms
ft(20001) = 1;                  % Define impulse at t = 2ms

fe = fx.'*ft; 

%% Processing 

% First time step is treated separately 
io(1) = accumFlux(0,iobar(1),T); 
id(1) = accumFlux(0,idbar(1),T); 
[csat(1), cs(1)] = saturate(ca(1), C, io(1), id(1)); 
[iobar(2), idbar(2)] = receptorCircuit(csat(1), io(1), id(1), kco, koc, kod, kdo, kcd, kdc);

ybar(:,1) = T*fe(:,1); 

for k = 2:length(t)-1
    
    % calculate boundary flux 
    ib(k) = iobar(k) + idbar(k);
    phib(:,k) = ca2_a.'*ib(:,k);
    
    % process inner model (state equation)
    ybar(:,k) = KeCe*A*ybar(:,k-1) + T*fe(:,k) - T*phib(:,k);
    
    % NT concentration at post synapse x = a
    ca(k) = c1_a*ybar(:,k);
    
    %%% Process Boundary Circuit %%% 
    
    % accumulate fluxes 
    io(k) = accumFlux(io(k-1),iobar(k),T); 
    id(k) = accumFlux(id(k-1),idbar(k),T);
    
    % saturate 
    [csat(k), cs(k)] = saturate(ca(k), C, io(k)/T, id(k)/T); 
    
    % number of receptors for next time step
    [iobar(k+1), idbar(k+1)] = receptorCircuit(csat(k), io(k), id(k), ... 
        kco, koc, kod, kdo, kcd, kdc);
end


%% Plot and comparison to PBS 
% addpath('./PBS-data')
% foo = readNPY('data_szenario2.npy');
% 
% occu = [0 cumsum(foo(:,1)).'];
% desen = [0 cumsum(foo(:,4)).'];
% open = occu - desen;
% 
% sample = round(length(t)/length(desen)); 
% 
% % Comparison
% figure(1); 
% plot(t_(1:sample:end)*1e3,desen,'or'); 
% grid on; hold on; plot(t_*1e3,id/T,'-b','LineWidth',1.5);
% xlabel('Time');
% ylabel('Number of Desensitized Receptors'); 
% legend('PBS', 'TFM');
% 
% figure(2); 
% plot(t_(1:sample:end)*1e3,open,'or'); 
% grid on; hold on; plot(t_*1e3,io/T,'-b','LineWidth',1.5);
% xlabel('Time');
% ylabel('Number of Open Receptors'); 
% legend('PBS', 'TFM');
% 
% figure(3); 
% plot(t_(1:sample:end)*1e3,occu,'or'); 
% grid on; hold on; plot(t_*1e3,io/T + id/T,'-b','LineWidth',1.5);
% xlabel('Time');
% ylabel('Number of Occupied Receptors'); 
% legend('PBS', 'TFM');

%% Plot stuff 

figure(1); plot(t_*1e3,io/T,'LineWidth',1.5);grid on;hold on; 
xlim([0 4])
title('Number of Receptors in OPEN state')
xlabel('Time in ms')
ylabel('Number of open Receptors');
legend('Setting 2.1', 'Setting 2.2', 'Setting 2.3', 'Setting 2.4')

figure(2); plot(t_*1e3,id/T,'LineWidth',1.5);grid on;hold on;
xlim([0 4])
title('Number of Receptors in DESENSITIZED state')
xlabel('Time in ms')
ylabel('Number of desensitized Receptors');
legend('Setting 2.1', 'Setting 2.2', 'Setting 2.3', 'Setting 2.4')

