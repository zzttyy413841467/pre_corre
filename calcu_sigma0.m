function sigma0=calcu_sigma0

constant_sim;
H0=121900;
V0=7630;
Hf=26500;
Vf=900;
r0=(Re+H0)/Re;
rf=(Re+Hf)/Re;
v0=V0/Vc;
vf=Vf/Vc;

e0=1/r0-v0^2/2;
ef=1/rf-vf^2/2;

espan=linspace(e0,ef,2000);
s0=68.5/57.3;
gamma0=-0.5*pi/180;
tau0=0;
y0=[s0;r0;gamma0;tau0];

sigma=zeros(2000,1);
qdot=zeros(2000,1);
sigma(1)=10/180*pi;
sigma(2)=15/180*pi;
t=1;
for i=1:1999
    sigma1=sigma(i);
    [ee,yy]=rk1(@dxde2,espan,y0,sigma1);
    ef=ee(end);
    rf=yy(end,2);
    rho=rho0.*exp(-(rf*Re-Re)./hs);
    v=sqrt(2*(1./rf-ef));
    qdot(i)=k_q*sqrt(rho)*v^3.15;
    qdot_delta=qdot(i)-q_max+10000;
    if abs(qdot_delta)<0.001
        break
    end
    
    if i>1
       sigma(i+1)=magnitude(sigma(i)-(sigma(i)-sigma(i-1))/(qdot(i)-qdot(i-1))*qdot_delta);   
    end
    t=t+1;
end
sigma0=180*sigma(t)/pi;

end