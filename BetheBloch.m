function [Stoppingpower]= BetheBloch(Z,A,I,rho,z,a,m,E,C,delta,VisBool)
%%%%%%%%%%%%%%%%%%%%
%Calculates the stoppingpower for a given atomicnumber z and kinetic energy in MeV/u of the projectile
%Z     charge of the target
%A     Atomicnumber of the target
%I     mean ionisation energy of the target
%rho   density of the target
%z     cahrge number of the projectile
%E     projectile energy in MeV/u
%C     shell correction
%delta density correction 
%m     projectile mass in atom
%
% refernece: https://doi.org/10.1063/1.369844
%%%%%%%%%%%%%%%%%%%%%


%% constants
c       = 299792458;               % Speed of light      m/s 
e       = 1.6021766208 * 10^(-19)  % electron charge     C
u       = 1.66053904   * 10^(-27)  % atomar mass         kg
epsilon = 8.854187817  * 10^(-12)  % epsilon nought      A*s/V/m
m_e     = 9.10938356   * 10^(-31)  % electron mass       kg
N_A     = 6.022140857  * 10^23     % avogadro constant   
alpha   = 1/137
h       = 6.626070040*10^(-34)/(2*pi)
r0      = e^2/(m_e*c^2);%0.52917721067 * 10^(-10) % bohr electron radius          m 
k       = 4 .* pi .* r0.^2 .* m_e * c.^2/e/1000/1000
%% Unit conversion

E=E*e*1000000 % energy in si unit

%% Bethe Bloch Formula

beta  = @(E,m) (sqrt(1-(m*c^2./(E+m*c^2)).^2))    % Formula for calculating beta from particle energie
gamma = @(E,m) (sqrt(1./(1-beta(E,m).^2)))
Tmax  = @(E,m) (2.*m_e.*c^2.*beta(E,m).^2 .* gamma(E,m).^2./((1+2.*gamma(E,m).*m_e./m+(m_e./m).^2)))    % Formula for retrieving Tmax for BB2





Stoppingpower =zeros(size(m,2)*size(Z,2),size(E,2));

for j=1:size(Z,2)
    for i = 1:size(m,2)
        k = 0.3071/Z(j)
        vBeta = beta(E*a(i),m(i));
        vL0 = (0.5 .*log(2.*m_e.*c.^2.*beta(E,m(i)).^2.*gamma(E,m(i)).^2.*Tmax(E,m(i))./I.^2)-beta(E,m(i)).^2-delta./2-C./Z.^2)     % Bethe-Bloch initial formula
        vLH = 1.5.*(E./e./1000).^(-0.4) + 45000./Z(j) .*(E./e./1000).^(-1.6);  
        vLL = 0.001*E/e/1000; 
        vL1 = vLL.*vLH./(vLH+vLL)                                                                                                   % Barkas-correction
        y   = z(i)*alpha./vBeta;
        vL2 =  -y.^2 .*(1.202-y.^2 .*(1.043-0.855 .* y.^2 + 0.343 .* y.^4))                                                         % Bloch-correction
        z_eff = z(i).*(1-exp(-125 .* beta(E*a(i),m(i))*z(i).^(-2/3)));                                                              % Netcliff z reduction for slow particles
        Stoppingpower(i*j,:) = (k .* Z(j)./vBeta.^2) .* z(i).^2 .* (vL0+vL1+vL2);
    end
end
E = E/e/1000000
if VisBool
    for j = 1:size(Z,2)
        figure%loglog(energie,dEdx.H.dEdx)
        strLegend =[]
        for i =1:size(m,2)
            loglog(E,Stoppingpower(i*j,:))
            hold on
            strLegend = [strLegend string(['Projectile with atomicmass ' num2str(m(i)/u) 'u and charge z=' num2str(z(i))])]
           
        end
        strtitle = ['Stoppingpower for target ' num2str(j)]
        legend(strLegend)
        title(strtitle)
        xlabel('Projectil Kinetic Energy [MeV/u]')
        ylabel('Stoppingpower [MeV/cm]')
        set(gca,'FontSize',14)
    
    end
end
end