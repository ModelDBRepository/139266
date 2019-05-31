y0 = [-67.7966    0.0661    0.9804    3.8280   20.0001    0  138.7929  143.9961    6.0  130.0]; % steady state

p.Tanoxia = 500*1e3;  %onset of anoxia
p.Ap = 0; p.Ad = .0;
p.Tcurr = [100,101]*1e3; % time between current is injected
p.Icurr = 1.6; % [uA/cm^2]
T = 2000*1e3; %1000*1e3; %[ms]

[Sol,p] = mainmodel(T,y0,p);

vars = fieldnames(p);
for i = 1:length(vars)
    assignin('base', vars{i}, p.(vars{i}));
end

%%
%teval = 0:1:200e3;             %solution can be resampled if desired
%y = zeros(10,length(teval));
%y = deval(teval,Sol);
%y = deval(Sol,10)

y = Sol.y';
t = Sol.x;

%%
%recalculate currents during simulation
Apump = ((Sol.x<Tanoxia)*(1-Ap)+Ap)'; 
alpha_m = 0.1 * (y(:,1)+30.0)./( 1.0 - exp(-0.1 * (y(:,1)+30.0)) );
beta_m = 4.0 * exp(-(y(:,1)+55.0)/18.0);
m_inf = alpha_m./(alpha_m + beta_m);
E_k = 26.64 * log((y(:,4)./y(:,7)));
E_na = 26.64 * log((y(:,8)./y(:,5)));
E_cl = -26.64 .* log((y(:,10)./y(:,9)));
Ina = g_na.*(m_inf.^3).*y(:,3).*(y(:,1)-E_na) + g_naL*(y(:,1)-E_na);
Ik = (g_k*y(:,2).^4 ).*(y(:,1)-E_k) + g_kL*(y(:,1)-E_k);
Icl = g_clL*(y(:,1)-26.64*log(y(9)./y(10)));
Ipump = Apump.*(rho./(1.0+exp((25.0-y(:,5))./3.0))).*(1./(1+exp(5.5-y(:,4))));

%%
%plot figures
y(y>1e10) = NaN;
figure(1); h1=axes; plot(t/1000,y(:,1)); xlabel('time (s)'); ylabel('Vm (mV)')
figure(2); plot(t/1000,[y(:,2),y(:,3)]); xlabel('time (s)'); title('opening n and h gates')
figure(3); h3=axes; plot(t/1000,[y(:,4),y(:,8),y(:,9),y(:,10)]); xlabel('time (s)'); ylabel('concentration mM');
figure(3); hold on; plot(t/1000,[y(:,7),y(:,5)],'--')
legend('[K]_e','[Na]_e','[Cl]_i','[Cl]_e','[K]_i','[Na]_i')
figure(4); h4=axes; plot(t/1e3,[E_k,E_na,y(:,1),E_cl]);  legend('E_k','E_{na}','V_m','E_{Cl}'); ylabel('V_m and Nernst potentials (mV)'); xlabel('time (s)')
figure(5); h5=axes; plot(t/1e3, [y(:,2),y(:,3)]); legend('n','h'); xlabel('time (s)'); title('K+ gate variables')

linkaxes([h1;h4;h3;h5],'x')