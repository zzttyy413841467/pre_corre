function [e,x1]=rk1(dxde2,espan,x0,sigma)
n=length(espan);
ef=espan(end);
e0=espan(1);
h=(ef-e0)/(n-1);
i=1;
x=x0;
x1=zeros(n,length(x0));
x1(1,:)=x;

for e=espan(1:end-1)
    
    result=trans(e,x);
    if (result==1)
        break;
    end
    k1=dxde2(e,x,sigma);
    k2=dxde2(e+h/2,x+h*k1/2,sigma);
    k3=dxde2(e+h/2,x+h*k2/2,sigma);
    k4=dxde2(e+h,x+h*k3,sigma);
    x=x+(h/6)*(k1+2*k2+2*k3+k4);
    i=i+1;
    x1(i,:)=x;
end
e=espan(1:i)';
x1(((i+1):n),:)=[];
end