function [Sol, p] = mainmodel(T,y0,p)
	%{
	y(1)=V, the membrane voltage
	y(2)=n, gating variable
	y(3)=h, gating variable
	y(4)=[K]_o, the extracellular potassium concentration
	y(5)=[Na]_i, the intracellular sodium concentration
    y(7)=[K]_i, the intracellular potassium concentration
    y(8)=[Na]_o, extracellular sodium concentration
    y(9)=[Cl]_i, intracellular chlorine concentration
    y(10)=[Cl]_o, extracellular chlorine concentration
	%}
    if nargin < 2       % set y0 if not given as input
        y0 = zeros(1, 10);

        %/* set initial conditions */
        y0(1)=-50;          %[mV]
        y0(2)=0.08553;      %[1]
        y0(3)=0.96859;      %[1]
        y0(4)=7.8;          %[mM]
        y0(5)=15.5;         %[mM]
        y0(6)= NaN;         %(not used in this simulation)
        y0(7)=140;          %[mM]
        y0(8)=144;          %[mM]
        y0(9)= 6;           %[mM]
        y0(10) = 130;       %[mM]
    end
    
	%/* set parameters */
    p.rcell = 7e-6;           % [m], radius of spherical cell
    p.F = 96485.3399;         % [C/mol], Faraday constant
    p.gamma = 1e-2*3/p.rcell/p.F; % [mM cm^2 /(uA s)] conversion from current to concentration change, gamma = A/(F*V) = 3/(rcell*F)
	p.tau = 1e3;           % conversion factor seconds -> ms
	p.beta = 2.0;             % ratio intra/extracellular volume;

	p.rho = 1.25/p.gamma;       % 1.25 mM/s / (mM cm^2 /(uA s)) = uA/cm^2 , pump current scaling
	p.glia = 200.0/3.0;       % mM/s, "pump rate" of [K+]e by glial cells
	p.epsilon = 4.0/3.0;      % [1/s] diffusion rate 
	p.kbath = 4.0;            % [mM], concentration K+ of "bath"

	p.Cm = 1.0;               % [uF / cm^2],  membrane capacitance
	p.g_na = 100.0;           % [mS / cm^2],  maximum gate conductances
	p.g_naL = 0.0175;         % [mS / cm^2],  leak conductance
	p.g_k = 40.0;             % [mS / cm^2]
	p.g_kL = 0.05;            % [mS / cm^2]
	p.g_clL = 0.05;           % [mS / cm^2]
	p.phi = 3.0;              % [1/ms],       gate time constant
    
    %Time in ms!
    
    %use odeprog.m to allow for interruption of calculation without loss of data
    %from www.mathworks.com/matlabcentral/fileexchange/9904-ode-progress-bar-and-interrupt  
    %tspan = [0,T]; options = odeset('MaxStep',1e3,'RelTol',1e-3,'OutputFcn',@odeprog,'Events',@odeabort);
    tspan = [0,T]; options = odeset('MaxStep',1e3,'RelTol',1e-3);
    funhan = (@(t,y)FullModel(t,y,p));
    Sol = ode23s(funhan,tspan,y0,options); %  Solve the differential equations using parameter set p
   
end

%---------------

function dydx = FullModel(t,y,p)

if t > p.Tanoxia;
    Apump = p.Ap; Adiff = p.Ad; Clconst = false;
else
    Apump = 1; Adiff = 1; Clconst = true;
end

if t > p.Tcurr(1) && t<p.Tcurr(2)   %inject current when specified
    Iapp = p.Icurr; %[uA/cm^2]
else
    Iapp = 0;
end

%x = time in ms

% Gates
	alpha_n = 0.01 * (y(1)+34.0)/( 1.0 - exp(-0.1 * (y(1)+34.0)) ); %[no units]
    beta_n = 0.125 * exp(-(y(1)+44.0)/80.0);
	alpha_m = 0.1 * (y(1)+30.0)/( 1.0 - exp(-0.1 * (y(1)+30.0)) );
	beta_m = 4.0 * exp(-(y(1)+55.0)/18.0);
	alpha_h = 0.07 * exp(-(y(1)+44.0)/20.0);
	beta_h = 1.0/( 1.0 + exp(-0.1 * (y(1)+14.0)) );
	m_inf = alpha_m/(alpha_m + beta_m);   
% Nernst potentials
    E_k = 26.64 * log(y(4)/y(7));       %[mV]
    E_na = 26.64 * log((y(8)/y(5)));    
    E_cl = 26.64*log(y(9)/y(10));
% Currents    
	Ina = p.g_na*(m_inf^3)*y(3)*(y(1)-E_na) + p.g_naL*(y(1)-E_na);  % [mS/cm^2 * mV = uA/cm^2]
	Ik = (p.g_k*y(2)^4)*(y(1)-E_k) + p.g_kL*(y(1)-E_k) ; 
    Icl = p.g_clL*(y(1)-E_cl);
	Ipump = Apump*(p.rho/(1.0+exp((25.0-y(5))/3.0)))*(1/(1+exp(5.5-y(4))));   % [mM/s]
	Iglia = Apump*(p.glia/(1.0+exp((18.0-y(4))/2.5)));                        % [mM/s]
	Idiffusion = Adiff*p.epsilon*(y(4)-p.kbath);                              % [mM/s]

	%{
	y(1)=V, the membrane voltage
	y(2)=n, gating variable
	y(3)=h, gating variable
	y(4)=[K]_o, the extracellular potassium concentration
	y(5)=[Na]_i, the intracellular sodium concentration
	y(6)=[Ca]_i, the intracellular calcium concentration
    y(7)=[K]_i, the intracellular potassium concentration
    %y(8)=[Na]_o, extracellular sodium concentration
    y(9)=[Cl]_i, intracellular chlorine concentration
    y(10)=[Cl]_o, extracellular chlorine concentration
	%}
    dydx = zeros(10,1);
	dydx(1) = (1/p.Cm)*(-Ina -Ik -Icl-0*Ipump+Iapp);
	dydx(2) = p.phi*(alpha_n*(1-y(2))-beta_n*y(2));
	dydx(3) = p.phi*(alpha_h*(1-y(3))-beta_h*y(3));
	dydx(4) = (1/p.tau)*(p.gamma*p.beta*Ik -2.0*p.beta*p.gamma*Ipump -Iglia -Idiffusion);
	dydx(5) = (1/p.tau)*(-p.gamma*Ina -3.0*p.gamma*Ipump);
	dydx(6) = 0;
    dydx(7) = -(1/p.tau)*(p.gamma*Ik -2.0*p.gamma*Ipump);
    dydx(8) = (1/p.tau)*(p.gamma*p.beta*Ina +3.0*p.beta*p.gamma*Ipump);
    if Clconst
        dydx(9) = 0;
        dydx(10)= 0;
    else
        dydx(9) = (1/p.tau)*(p.gamma*Icl);
        dydx(10)= -dydx(9)*p.beta;
    end
end