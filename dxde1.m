function xdot=dxde1(e,x)
constant_sim;

alpha0=45/180*pi;
sigmaf=evalin('base','sigmaf');
sigma0=evalin('base','sigma0');
ef=evalin('base','ef');
e0=evalin('base','e0');
sigma=limit(sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0),e,x);
% % sigma=sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0)

s=x(1);
r=x(2);
gamma=x(3);
V=sqrt(2*(1./r-e));

ma=V*Vc/340;
alpha=alpha0.*(ma>=10)+((45-0.612*(ma-10).^2)/180*pi).*(ma<10);

Cl = cl0+cl1*alpha+cl2.*alpha.^2;
Cd = cd0+cd1*Cl+cd2.*Cl.^2;

h=r*Re-Re;
rho=rho0.*exp(-h./hs);
q=1/2*rho.*V.^2;
L=q.*Cl.*S/m*Re;
D=q.*Cd.*S/m*Re;
% g=mu./(Re+h).^2;

rdot=sin(gamma)./D;

gammadot=1./(D.*V.^2).*(L.*cos(sigma)+(V.^2./r-1./r.^2).*cos(gamma));
sdot=-cos(gamma)./(r.*D);
taudot=1./(D.*V);

theta=x(4);
phi=x(5);
psi=x(6);
thetadot=cos(gamma).*sin(psi)./(r.*D.*cos(phi));
phidot=cos(gamma).*cos(psi)./(r.*D);
psidot=1./(D.*V).*(L.*sin(sigma)./V./cos(gamma)+V./r.*cos(gamma).*sin(psi).*tan(phi));

xdot=[sdot;rdot;gammadot;taudot;thetadot;phidot;psidot];

end