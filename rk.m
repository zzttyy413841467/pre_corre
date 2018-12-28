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
sigma(1)=evalin('base','sig0')/180*pi;
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
    xx0=[0,x(1),x(4)];
    sigma_x=limit(sigma0+(sigmaf-sigma0)./(ef-e0).*(e-e0),e,xx0);
    assignin('base','sigma_x',sigma_x);
    [ee1,xx1]=ode45(@dyde,[e,ef],x(1:5));
    sf(i)=acos(cos(xx1(end,3))*cos(phif)*cos(xx1(end,2)-thetaf)+sin(xx1(end,3))*sin(phif));
    
    
    if i>=2
       sigma(i)=sigma_x*sign_decide(sign(sigma(i-1)),e,x); 
    end
    k1=dxde1(e,x,sigma(i));
    k2=dxde1(e+h/2,x+h*k1/2,sigma(i));
    k3=dxde1(e+h/2,x+h*k2/2,sigma(i));
    k4=dxde1(e+h,x+h*k3,sigma(i));
    x=x+(h/6)*(k1+2*k2+2*k3+k4);
    i=i+1;
    x1(i,:)=x;
end
e=espan';
end