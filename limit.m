function sigma=limit(sigma1,e,x)

constant_sim;

sigma0=20*pi/180;

r=x(2);
v=sqrt(2*(1/r-e));
gamma=x(3);
alpha0=45/180*pi;
ma=v*Vc/340;
alpha=alpha0.*(ma>=10)+((45-0.612*(ma-10).^2)/180*pi).*(ma<10);

Cl = cl0+cl1*alpha+cl2.*alpha.^2;
Cd = cd0+cd1*Cl+cd2.*Cl.^2;

Kl=Re*Cl*S/2/m;

h=r*Re-Re;
rho=rho0.*exp(-h./hs);
q=1/2*rho.*v.^2;
L=q.*Cl.*S/m*Re;
D=q.*Cd.*S/m*Re;
drhodr=rho*(-Re/hs);
drdv_q=-6.3*q_max^2/(k_q^2*drhodr*v^7.3);
drdv=v.*sin(gamma)./(-D-sin(gamma)./r.^2);
eps0=0*(drdv-drdv_q);

sigmaqdot=acos(k_q.^2.*(1-v.^2+eps0).*v.^4.3./(Kl*q_max^2));
sigmaq=acos(Vc^2*(1-v.^2)./(2*Kl*1.5e4));
sigman1=acos((1-v.^2)./2.5.*(cos(alpha)+Cd./Cl.*sin(alpha)));

if sigma1>min([sigmaqdot sigmaq sigman1])
    sigma1=min([sigmaqdot sigmaq sigman1]);
end

if sigma1<sigma0
    sigma1=sigma0;
end
sigma=sigma1;


end