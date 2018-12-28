function ydot=dyde(e,y)

constant_sim;

alpha0=45/180*pi;
sigma0=evalin('base','sigma0');
sigmaf=evalin('base','sigmaf');
e0=evalin('base','e0');
ef=evalin('base','ef');
sigma=limit(sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0),e,y);
% sigma=sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0)
r=y(1);
gamma=y(2);
V=sqrt(2*(1./r-e));

ma=V*Vc/340;
alpha=alpha0.*(ma>=15)+((45-0.21*(ma-15).^2)/180*pi).*(ma<15);

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
ydot=[rdot;gammadot];


end