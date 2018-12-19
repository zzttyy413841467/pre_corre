function xdot=dxde(e,x)

constant_sim;

alpha0=45/180*pi;
% sigma0=evalin('base','sigma0');
% sigmaf=evalin('base','sigmaf');
% e0=evalin('base','e0');
% ef=evalin('base','ef');
% sigma=sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0);

r=x(1);
theta=x(2);
phi=x(3);
gamma=x(4);
psi=x(5);
v=sqrt(2*(1./r-e));

ma=v*Vc/340;
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
thetadot=cos(gamma).*sin(psi)./(D.*r.*cos(phi));
phidot=cos(gamma).*cos(psi)./(D.*r);
gammadot=1./(D.*v.^2).*(L.*cos(sigma)+(v.^2./r-1./r.^2).*cos(gamma)+2*omega.*v.*cos(phi).*sin(psi)+omega^2.*r.*cos(phi).*(cos(gamma).*cos(phi)+sin(gamma).*cos(psi).*sin(phi)));
psidot=1./(D.*v.^2).*(L.*sin(sigma)./cos(gamma)+(v.^2./r).*cos(gamma).*sin(psi).*tan(phi)-2*omega.*v.*(tan(gamma).*cos(phi).*cos(psi)-sin(phi))+omega^2.*r.*cos(phi).*sin(phi).*sin(psi)./cos(gamma));
% sdot=-cos(gamma)./(r.*D);
taudot=1./(D.*v);
xdot=[rdot;thetadot;phidot;gammadot;psidot];


end