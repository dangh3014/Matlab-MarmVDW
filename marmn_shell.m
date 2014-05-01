% Nondimensional Marmottant bubble solver
% To dimensionalize: R = R*R0; V = V*R0*w

% tic
% matlabpool                    % Parallelize
% cd /data/bubble/daking3/simulations/

clear all


%%%%% Bubble Properties %%%%%
R0_vect = 2e-6;                 % [m]
R_pert = 0;                     % Initial perturbation (i.e. 0.9, 1.001) (Pulse will be ignored); if 1, will run as pert.
V0 = 0;                         % [m/s] Initial velocity
gas = 'vdw';                    % 'vdw' or 'ideal'
bubble_eq = 'free';             % 'marm', 'dejong', or 'free'

%%%%% Shell Properties %%%%%
% Note: for de Jong model, Chi = Sp/2, KappaSh = Sf/(16*pi) (<- Note this is NOT Sf/(12*pi) as in paper - see Overvelde dissertation)
KappaSh = 2.4e-9;               % Optison (Chatterjee): 7.65e-8; Definity (Goertz): 0.398e-9; Definity2 (Kimmel): 2.4e-9; Free = 0;
Chi = 0.38;                     % Optison (Chatterjee): 0.932; Definity (Goertz): 0.855; Definity2 (Kimmel): 0.38; Free = 0;

Rbuck = 0.99;                   % Normalized to R0 = 1;
SigmaR0 = Chi*((1/Rbuck)^2-1);

SigmaL = 0.073;                 % [N/m] Surface tension of liquid
Rrupt = Rbuck*(1+SigmaL/Chi)^(1/2);
Rbreak = Rrupt;                 % 1.5 1.2 1.08 Rrupt 1.3

% SigmaBreak = 50;                % Unused
% Rbreak = Rbuck*(1+SigmaBreak/Chi)^(1/2); % Theoretical


%%%%% Pulse Properties %%%%%
pulse_type = 'windowed';        % 'sin', 'windowed', or 'kzk'
t=(0:0.001:5)*10^-6;            % [s] Simulation time

PRP_vect = 0.1e6;                 % [Pa]
F_vect = 2.8e6;                   % [Hz]
Cycles = 3;                     % # of cycles
         
        
% Fluid and Gas Parameter Values
Rho = 998;       % [kg/m] Density of liquid
P0 = 101325;     % [Pa] Hydrostatic Pressure 
SigmaL = 0.073;  % [N/m] Surface tension of liquid
C = 1481;        % [m/s] Velocity of Sound in liquid
Mu = .001;       % [Pa*s] (Dynamic) Viscosity of surrounding liquid (Wu suspension 1.28e-3)
                 % Nu = Kinematic Viscosity = mu/rho

KappaG = 1.06;   % Polytropic Exponent, Definity C3F8)
                 % Definity 1.06 (Goerze) 
                 % 1.095 SF6
                 % 1.4 Air 
                 % Optison 1.09
                 % BEM 1.25

if strcmp(pulse_type,'kzk')

    load /data/bubble/daking3/simulations/kzk_pulses/kzk_5MHz_3cycle

    R=zeros(size(p,2),size(p,1),length(R0_vect)); % 49, 99
    V=zeros(size(p,2),size(p,1),length(R0_vect)); % 49, 99

else               
                 
    R=zeros(length(t),size(PRP_vect(:),1),size(R0_vect(:),1),size(F_vect(:),1));
    V=zeros(length(t),size(PRP_vect(:),1),size(R0_vect(:),1),size(F_vect(:),1));

end

for jj=1:size(PRP_vect(:),1)
    
for kk=1:size(F_vect(:),1)
       
    F = F_vect(kk);

    w=2*pi*F;
    T=t*w;
    
    if R_pert==0               % Simulate Pulses

    R_pert_ini = 1;
    Period = 2*pi;
    T0 = Period*Cycles;   

switch pulse_type
    
    case 'sin'
    
    pulse_window = zeros(1,length(T));
    for m=1:length(T)
        if T(m)>0 && T(m)<T0
            pulse_window(m) = 1;
        else
            pulse_window(m) = 0;
        end
    end
        
    case 'windowed' % To simulate transducer ring up and ring down
        
    pulse_window = zeros(1,length(T));
    for m=1:length(T);
        if T(m)>0 && T(m)<T0
            pulse_window(m) = (1-exp(-T(m)/Period))/(exp(-T0/Period)-1);
        elseif T(m)>=T0 && T(m)<5*T0
            pulse_window(m) = -exp(1-T(m)/Period)/exp(1-T0/Period);
        else
            pulse_window(m) = 0;
        end
    end

end

if strcmp(pulse_type,'kzk')
    P=p(jj,:,2)*10^6/P0;
else
    PRP = PRP_vect(jj);
    A = PRP/P0;                 % Nondimensional
    P = -A*sin(w*t).*pulse_window;    
end

    else                        % Perturbation Only

R_pert_ini = R_pert;
P=zeros(1,length(t)); 
    
    end

for ll=1:length(R0_vect)

% % Growth at breaking point
% s = 1.4;            
% temp = find(Rnd(:,j,19)>=s,1,'first');
% if isempty(temp);
%     continue
% else
% bp(j,k) = temp;
% Rini = [Rnd(bp(j,k),j,19) Vnd(bp(j,k),j,19)];
% tspan = [T(bp(j,k)) T(end)];

R0=R0_vect(ll); 
Rini = [R_pert_ini V0]';
tspan = [T(1) T(end)];      % To cut out 1st us of KZK: [T(671) T(end)] (deval from 671:end)

options = [];               %['RelTol',1e-6,'AbsTol',1e-9]; %[]odeset('OutputFcn',@odeplot); %OutputFcn', @odeplot, @hdot ,'MaxStep',.1
    
% Useful ode solvers: ode 45, 15s, 23s
Rn = ode23s(@marmottantn_vdw,tspan,Rini,options,T,P,w,R0,Rbuck,Rrupt,Rbreak,KappaSh,Chi,SigmaR0,Rho,P0,SigmaL,C,Mu,KappaG,bubble_eq,gas);

if max(Rn.x)>=max(T)
R(:,jj,kk,ll) = deval(T,Rn,1);
V(:,jj,kk,ll) = deval(T,Rn,2);
else        % if solver breaks before end
tempt=find(T<max(Rn.x),1,'last');
R(1:tempt,jj,kk,ll) = deval(T(1:tempt),Rn,1);
V(1:tempt,jj,kk,ll) = deval(T(1:tempt),Rn,2);
end



end
end
end

properties = struct('R0_vect',R0_vect,'R_pert',R_pert,'V0',V0,'gas',gas,'bubble_eq',bubble_eq,...
    'KappaSh',KappaSh,'Chi',Chi,'Rbuck',Rbuck,'SigmaR0',SigmaR0,'SigmaL',SigmaL,'Rrupt',Rrupt,...
    'Rbreak',Rbreak,'pulse_type',pulse_type,'PRP_vect',PRP_vect,'F_vect',F_vect,'Rho',Rho,...
    'P0',P0,'C',C,'Mu',Mu,'KappaG',KappaG,'ode_solver','ode23s');

% Typical vars to save: R V R0_vect F t prp properties
% save -v7.3 marm_sim R V R0_vect F t prp properties



% matlabpool close
% toc
