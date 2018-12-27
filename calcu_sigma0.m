function calcu_sigma0

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

espan=linspace(e0,ef,5000);
gamma0=-0.5*pi/180;
y0=[r0;gamma0];

sigma=zeros(2000,1);
qdot=zeros(2000,1);
sigma(1)=30/180*pi;
sigma(2)=35/180*pi;
t=1;
for i=1:1999
    assignin('base','sig0_a',sigma(i)*180/pi);
    options=odeset('events',@stop_conditions_1);
    [ee,yy]=ode45(@dxde2,[e0 ef],y0,options);
    ef=ee(end);
    rf=yy(end,1);
    rho=rho0.*exp(-(rf*Re-Re)./hs);
    v=sqrt(2*(1./rf-ef));
    qdot(i)=k_q*sqrt(rho)*v^3.15;
    qdot_delta=qdot(i)-q_max;
    
    if abs(qdot_delta)<0.001
        break
    end
    
    if i>1
       sigma(i+1)=magnitude(sigma(i)-(sigma(i)-sigma(i-1))/(qdot(i)-qdot(i-1))*qdot_delta);   
    end
    t=t+1;
end

end