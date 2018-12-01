function [e,x1,sigma]=rk(dxde1,espan,x0)
n=length(espan);
ef=espan(end);
e0=espan(1);
h=(ef-e0)/(n-1);
i=1;
x=x0;
x1=zeros(n,length(x0));
sigmaf=evalin('base','sigmaf');
sigma_0=zeros(n,1);
sigma=zeros(n-1,1);
x1(1,:)=x;
sigma_0(1)=evalin('base','sigma0');
sigma0=sigma_0(1);
sf=zeros(n,1);
thetaf=evalin('base','thetaf');
phif=evalin('base','phif');
for e=espan(1:end-1)
    if i==2
        sigma0=sigma0+10/57.3;
        sigma_0(i)=sigma0;
        assignin('base','sigma0',sigma0);
    end
    if i>2
        sigma0=magnitude(sigma_0(i-1)-(sigma_0(i-1)-sigma_0(i-2))/(sf(i-1)-sf(i-2))*sf(i-1));
        sigma_0(i)=sigma0;
        assignin('base','sigma0',sigma0);
    end
    
    xx0=[acos(cos(x(3))*cos(phif)*cos(x(2)-thetaf)+sin(x(3))*sin(phif)),x(1),x(4)];
    [ee1,xx1]=ode45(@dyde,linspace(e,ef,200),xx0);
    sf(i)=xx1(end,1);
    sigma_x=limit(sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0),e,xx0);
    
    sigma(i)=sigma_x;
    k1=dxde1(e,x);
    k2=dxde1(e+h/2,x+h*k1/2);
    k3=dxde1(e+h/2,x+h*k2/2);
    k4=dxde1(e+h,x+h*k3);
    x=x+(h/6)*(k1+2*k2+2*k3+k4);
    i=i+1;
    x1(i,:)=x;
end
e=espan';
end