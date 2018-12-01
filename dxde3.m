function xdot=dxde3(e,x)
constant_sim;

sigma=evalin('base','sig0')/180*pi;
alpha0=45/180*pi;

r=x(1);
gamma=x(4);
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
taudot=1./(D.*V);
theta=x(2);
phi=x(3);
psi=x(5);
thetadot=cos(gamma).*sin(psi)./(r.*D.*cos(phi));
phidot=cos(gamma).*cos(psi)./(r.*D);
psidot=1./(D.*V).*(L.*sin(sigma)./V./cos(gamma)+V./r.*cos(gamma).*sin(psi).*tan(phi));

xdot=[rdot;thetadot;phidot;gammadot;psidot;taudot];
end