% Full non-dimensional Marmottant
% 8/3/09

function [Rn] = marmottantn(t,R,T,P,w,R0,Rbuck,Rrupt,Rbreak,KappaSh,Chi,SigmaBreak,SigmaR0)

% Parameter Values
Rho = 1000;%998;       %[kg/m] Density of liquid
P0 = 101325; %1.07e8;    %[Pa] Hydrostatic Pressure 
SigmaL = 0.073;  %[N/m] Surface tension of liquid
C = 1480;        %[m/s] Velocity of Sound in liquid
Mu = .001;       %[Pa*s] (Dynamic) Viscosity of surrounding liquid (Wu suspension 1.28e-3)
                 %Nu = Kinematic Viscosity = mu/rho

KappaG = 1.25;   %Polytropic Exponent, Definity C3F8)
                 % Definity 1.06 (Goerze) 
                 % 1.095 SF6
                 % 1.4 Air 
                 % Optison 1.09
                 % BEM 1.25

% Pulse Interpolation
Pac = interp1(T,P,t);  %Interpolate pulse  


% Shell State (not allowing breakup or rupture - see marmottant1b,2)
if R(1) < Rbuck
   SigmaR = 0;

elseif R(1) >= Rbuck && R(1) < Rbreak %Rrupt
   SigmaR = Chi*((R(1)/Rbuck)^2-1);
   
elseif R(1) >= Rbreak %Rrupt 
   SigmaR = SigmaL;
   
else
   SigmaR = SigmaL;

end

% Free - uncomment this section
KappaSh = 0;
SigmaR = 0; % SigmaL;
SigmaR0 = 0; % SigmaL; 

% % De Jong-like - uncomment this section
% SigmaR = SigmaR0+2*Chi*((R(1)/Rbuck)-1);

%Solve Paired Equations

%Note: R(1) = Radius, R(2) = dR/dt
%Normalized radius Rn = R/R0 -> R = R0*Rn, etc
Rn(1) = R(2);
Rn(2) = (1/(Rho*w^2*R0^2*R(1)))*...
    ((-3/2*Rho*(R0*w*R(2))^2)+...
    (P0+2*SigmaR0/R0)*(1/R(1))^(3*KappaG)...
    -4*KappaSh*w*R(2)/(R0*R(1)^2)...
-2*SigmaR/(R0*R(1))... 
-(P0+P0*Pac));

%    *(1-3*KappaG/C*R0*w*R(2))...
%      -4*Mu*w*R(2)/R(1)...

Rn = Rn(:); %Column vectorize

end
