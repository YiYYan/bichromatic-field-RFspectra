clear all;
clc;
global delta kappa Ax Davg r LL sz sp sn w1 w2 w0
kappa=1.0;
Davg=5.0;
delta=0.1*kappa;
Ax=5*kappa;w0=250.0;w1=w0-0.1;w2=w0+0.1;
r=1.0;
sp=[0.0 1.0;0.0 0.0];
sz=[1.0 0.0;0.0 -1.0];
sn=sp.';
LL=0.5*kappa*(kron(eye(2),sp*sn)+kron(sp*sn,eye(2))-2*kron(sn,sn));
y0=[1;0;0;0];
tspan=linspace(0,300,3001);
nt=length(tspan);
% options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[tt0,yy0]=ode45(@dydtfun,tspan,y0);
% figure(2);plot(tt0,yy0(:,2)),hold on
gxxfun=diag(yy0(:,1));
for ii=1:nt-1
    yyt1=reshape(sn*reshape(yy0(ii,:),2,2),[],1);
    if ii==nt-1
        tspan2=tspan(ii:end);
        [tt2,yy2]=ode45(@(t,y) dydtfun(t,y),tspan2,yyt1);
        gxxfun(ii,ii:end)=yy2([1,end],2).';
    else
        tspan2=tspan(ii:end);
        [tt2,yy2]=ode45(@(t,y) dydtfun(t,y),tspan2,yyt1);
        gxxfun(ii,ii:end)=yy2(:,2).';
    end
end
[T1,T2]=meshgrid(tspan,tspan);
TF=tspan(end);
Gamma=0.0*kappa;
gxxfun=gxxfun+tril(gxxfun',-1);
BB=3001;
mm=length(BB);
wxl=linspace(-15,15,500);
Nw=length(wxl);
swt=zeros(mm,Nw);
swtc=zeros(mm,Nw);
for nn=1:mm
    AA=BB(nn);
    TT1=T1(1:AA,1:AA);
    TT2=T2(1:AA,1:AA);
    GXfun=gxxfun(1:AA,1:AA);
    for ii=1:Nw
        swt(nn,ii)=trapz(tspan(1:AA),trapz(tspan(1:AA),GXfun.*exp(-1i*wxl(ii)*(TT1-TT2)).*exp(-0.5*Gamma*(TF-TT1)).*exp(-0.5*Gamma*(TF-TT2)),2));
        swtc(nn,ii)=trapz(tspan(1:AA),(yy0(1:AA,2).').*exp(-1i*wxl(ii)*tspan(1:AA)));
    end
end
% plot(tspan(BB),swt)
% swtr=jwfun(wxl).*swt;
swtr=real(swt);
plot(wxl,swtr-abs(swtc).^2),hold on
% waterfall(swtr)