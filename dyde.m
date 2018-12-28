function ydot=dyde(e,y)

constant_sim;

alpha0=45/180*pi;
sigma0=evalin('base','sigma0');
sigmaf=evalin('base','sigmaf');
e0=evalin('base','e0');
ef=evalin('base','ef');
sigma=limit(sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0),e,[0,y(1),y(4)]);
sigma_x=evalin('base','sigma_x');
sigma=sign_decide(sign(sigma_x),e,y)*sigma;
assignin('base','sigma_x',sigma);
% sigma=sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0)
r=y(1);
gamma=y(4);
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

theta=y(2);
phi=y(3);
psi=y(5);
thetadot=cos(gamma).*sin(psi)./(r.*D.*cos(phi));
phidot=cos(gamma).*cos(psi)./(r.*D);
psidot=1./(D.*V).*(L.*sin(sigma)./V./cos(gamma)+V./r.*cos(gamma).*sin(psi).*tan(phi));


rdot=sin(gamma)./D;
gammadot=1./(D.*V.^2).*(L.*cos(sigma)+(V.^2./r-1./r.^2).*cos(gamma));
ydot=[rdot;thetadot;phidot;gammadot;psidot];


end