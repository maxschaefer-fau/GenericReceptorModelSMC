% General Receptor based on [1]
%
% M. Schäfer, S. Lotter, M. T. Barros, 2022 
% Further details can be found in: 
% 
%
% General 3-state receptor model (Open, Closed, Desensitized) with
% possible saturation from O->C and O->D. The model can be specialized to
% any 3-state receptor by inserting the desired transition rates. 
%
%
% References: 
%   [1] S. Lotter, M. T. Barros, R. Schober, M. Schaefer, "Signal Reception 
%       With Generic Three-State Receptors in Synaptic MC", submitted to
%       IEEE Global Commun. Conf. (GLOBECOM 2021), Arxiv: 
%   
%   [2] ﻿S. Lotter, J. Zeitler, M. Schaefer, R. Schober, "Saturating receiver 
%       and receptor competition in synaptic DMC: Deterministic and 
%       statistical signal models", IEEE Trans. Nanobiosci., 20(4), 464–479 
%       2021, https://doi.org/10.1109/TNB.2021.3092279.
%
%


%% Paramter Definition (with physical Dimensions taken from [2]) 

T_ = 1e-8;       % simulatiton time step in s 
t_end_ = 6e-3;   % simulation duration in s 
t_ = 0:T_:t_end_;  % time vector 

Na = 6.022e23; % Avogadro constant in mol^-1

a_ = 2*1e-2*1e-6; % channel width in x-direction in m
yz_ = 0.15;  % channel width in y- and z-direction in mum (1e-6 m) 

D_ = 3.3e-4*1e-12/1e-6; % Diffusion coefficient in m^2/s
keCe_ = 1e-3/1e-6; % Enzymatic degradation rate in 1/s

%% Receptor Definition (with physical Dimensions)
recepDensity_ = 1e3/1e-12;  % Receptor density in 1/m^2 TODO 
C = 203;                   % number of receptors 
r_rec_ = 2.3e-3;            % receptor radius in mum
rho_ = (C*pi*r_rec_^2)/yz_^2; % dimensionless receptor coverage


% Default setting for kinetic parameters ( see [1, Section IV-A]) 
% AMPAR
% C <---> O 
kco_ = 0.995*1.02e-4*rho_;   % in m/s
koc_ = 8.5e-3/1e-6;          % in 1/s
%O <---> D 
kod_ = 8.5e-3/1e-6;          % in 1/s
kdo_ = 0;                    % in 1/s
%C <---> D 
kcd_ = 0;                    % in m/s
kdc_ = 8.5e-3/1e-6;            % in 1/s

% NMDAR 
% C <---> O 
% kco_ = 0;   % in m/s
% koc_ = 8.5e-3/1e-6;          % in 1/s
% % O <---> D 
% kod_ = 0;          % in 1/s
% kdo_ = 8.5e-3/1e-6;                    % in 1/s
% % C <---> D 
% kcd_ = 0.995*1.02e-4*rho_;                    % in m/s
% kdc_ = 8.5e-3/1e-6;            % in 1/s

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
koc = koc_*a_^2/D_; 
kod = kod_*a_^2/D_; 
kdo = kdo_*a_^2/D_; 
kcd = kcd_*a_/D_; 
kdc = kdc_*a_^2/D_; 

%% Rate variations 
% The dimensionless rates kco-kdc can be varied to study the receptor
% behavior. 

% Vector of rate variation coefficients 
sc = [1   1  1  1      1   1]; 

[kco, koc, kod, kdo, kcd, kdc] = rateVariation(kco, koc, kod, kdo, ...
     kcd, kdc, sc);


% Define release pattern for NTs (single, or multiple)
exci = 'multiple'; 


%% Set up transfer function model for the inner synaptic model 

% see below eq. (26) in [2]
Mu = 10;           % number of eigenvalues 
mu = 0:Mu-1;        % vector to count eigenvalues 
gmu = mu*pi/a;      % wave numbers gamma 
smu = -D*gmu.^2;    % eigenvalues of the inner synaptic model 

% Eigenfunctions, see eq. (36) in [2]
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

t0 = 0; 
t1 = 1e-3; 
t2 = 2e-3;

if(isequal(exci,'single'))
    ft(1) = 1;                  % Define impulse at t = 0
else 
    ft(1) = 1;                        % Define impulse at t = 0
    ft(t1/T_+1) = 1;                  % Define impulse at t = t1
    ft(t2/T_+1) = 1;                  % Define impulse at t = t2
end


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
    
    % accumulate fluxes, see [1, Eq. (17)]
    io(k) = accumFlux(io(k-1),iobar(k),T); 
    id(k) = accumFlux(id(k-1),idbar(k),T);
    
    % saturate [1, Eq. (20)]
    [csat(k), cs(k)] = saturate(ca(k), C, io(k)/T, id(k)/T); 
    
    % number of receptors for next time step, see [1, Eqs. (18), (19)]
    [iobar(k+1), idbar(k+1)] = receptorCircuit(csat(k), io(k), id(k), ... 
        kco, koc, kod, kdo, kcd, kdc);
end


%% Plot stuff (Dimensionless Time and dimensional number of receptors)

figure(1); plot(t,io/T,'LineWidth',1.5);grid on;hold on; 
title('Number of Receptors in OPEN state')
xlabel('Dimensionless time')
ylabel('Number of open Receptors');

ax = get(gca);
ax.XAxis.Exponent = 3;

figure(2); plot(t,id/T,'LineWidth',1.5);grid on;hold on;
title('Number of Receptors in DESENSITIZED state')
xlabel('Dimensionless time')
ylabel('Number of desensitized Receptors');

ax = get(gca);
ax.XAxis.Exponent = 3;


