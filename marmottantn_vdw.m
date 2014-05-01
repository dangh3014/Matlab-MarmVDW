% Full non-dimensional Marmottant

function [Rn] = marmottantn_vdw(t,R,T,P,w,R0,Rbuck,Rrupt,Rbreak,KappaSh,Chi,SigmaR0,Rho,P0,SigmaL,C,Mu,KappaG,bubble_eq,gas)

% % Parameter Values
% Rho = 1000;%998;       %[kg/m] Density of liquid
% P0 = 101325;     %[Pa] Hydrostatic Pressure 
% SigmaL = 0.073;  %[N/m] Surface tension of liquid
% C = 1480;        %[m/s] Velocity of Sound in liquid
% Mu = .001;       %[Pa*s] (Dynamic) Viscosity of surrounding liquid (Wu suspension 1.28e-3)
%                  %Nu = Kinematic Viscosity = mu/rho
% 
% KappaG = 1.06;   % Polytropic Exponent, Definity C3F8)
%                  % Definity 1.06 (C3F8 Goerze) 
%                  % 1.095 SF6
%                  % 1.4 Air 
%                  % Optison 1.09

% Pulse Interpolation
Pac = interp1(T,P,t);  %Interpolate pulse  


% Shell State (not allowing breakup or rupture - see marmottant1b,2)
switch bubble_eq
    case 'marm'

    if R(1) < Rbuck
        SigmaR = 0;

    elseif R(1) >= Rbuck && R(1) < Rbreak   %Rrupt
        SigmaR = Chi*((R(1)/Rbuck)^2-1);

    elseif R(1) >= Rrupt                    %Rbreak
        SigmaR = SigmaL;

    else
        SigmaR = SigmaL;

    end
    
    case 'dejong'
    % De Jong-like - uncomment this section
    SigmaR0 = SigmaL;
    KappaSh = KappaSh*16*pi; % Not 12*pi!
    SigmaR = SigmaL+2*Chi*(R(1)*R0)*(1/R0-1/(R(1)*R0));

    case 'free'
    % Free - uncomment this section
    KappaSh = 0;
    SigmaR = SigmaL;
    SigmaR0 = SigmaR;
    
end


%Solve Paired Equations

switch gas
    
    case 'vdw'
% van der Waals hard core radius (C3F8)
h=R0/5.6; % R0/5.6; R0/5.644; R0/5.612; Air: R0/8.54 (Loefstedt)

%Note: R(1) = Radius, R(2) = dR/dt
%Normalized radius Rn = R/R0 -> R = R0*Rn, etc
Rn(1) = R(2);
Rn(2) = (1/(Rho*w^2*R0^2*R(1)))*...
    ((-3/2*Rho*(R0*w*R(2))^2)+...
    (P0+2*SigmaR0/R0)*(((R(1)*R0)^3-h^3)/(R0^3-h^3))^(-KappaG)*...
    (1-3*KappaG/C*(R0*w*R(2))*(R(1)*R0)^3/((R(1)*R0)^3-h^3))...
    -2*SigmaR/(R0*R(1))...
    -4*Mu*w*R(2)/R(1)...
    -4*KappaSh*w*R(2)/(R0*R(1)^2)...
    -(P0+P0*Pac));

    case 'ideal'
        
Rn(1) = R(2);
Rn(2) = (1/(Rho*w^2*R0^2*R(1)))*...
    ((-3/2*Rho*(R0*w*R(2))^2)+...
    (P0+2*SigmaR0/R0)*(1/R(1))^(3*KappaG)...
    *(1-3*KappaG/C*R0*w*R(2))...
    -2*SigmaR/(R0*R(1))... 
    -4*Mu*w*R(2)/R(1)...
    -4*KappaSh*w*R(2)/(R0*R(1)^2)...
-(P0+P0*Pac));
        
end

Rn = Rn(:); %Column vectorize

end
