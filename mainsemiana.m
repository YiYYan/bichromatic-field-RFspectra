clc;clear;
kx=2500;% for the first mode
ky=0;% for the second mode, if consider single mode drive, set ky=0; 
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
ui=1:2:Nh-1;
di=2:2:Nh;
[Imx,Dmx,Cmx]=genelemat(kx);
sx=[1 0;0 -1];
sz=[0 1;1 0];
sy=-[0 -1i;1i 0];
sp=(sx+1i*sy)/2;
sn=(sx-1i*sy)/2;
kappa=1.0;
Davg=1.0;
delta=0.001*kappa;
A=10.0*kappa;r=1;
H0=0.5*Davg*sz;
H1p=0.5*A*(r*sn+sp);
H1n=0.5*A*(sn+r*sp);
Hf=kron(Imx,H0)+kron(Dmx,delta*eye(2))+kron(Cmx,H1n)+kron(Cmx.',H1p);
[u,v]=eig(Hf);
rc=find(diag(v)<0.5*delta&diag(v)>=-0.5*delta);
er=diag(v(rc,rc));
ba=u(:,rc);
Der=er(1)-er(2);
f=@(a,b,k) sum(a(1-k*(k<0):Nk-k*(k>0)).*b(1+k*(k>0):Nk+k*(k<0)));
xppl=zeros(1,Nk);
xpnl=zeros(1,Nk);
xnpl=zeros(1,Nk);
for k=kxv
    xppl(k+kx+1)=f(ba(ui,1),ba(ui,1),k)*0.5-f(ba(di,1),ba(di,1),k)*0.5-f(ba(ui,1),ba(di,1),k)*0.5+0.5*f(ba(di,1),ba(ui,1),k);%f(ba(ui,1),ba(di,1),k);
    xpnl(k+kx+1)=f(ba(ui,1),ba(ui,2),k)*0.5-f(ba(di,1),ba(di,2),k)*0.5-f(ba(ui,1),ba(di,2),k)*0.5+0.5*f(ba(di,1),ba(ui,2),k);%f(ba(ui,1),ba(di,2),k);
    xnpl(k+kx+1)=f(ba(ui,2),ba(ui,1),k)*0.5-f(ba(di,2),ba(di,1),k)*0.5-f(ba(ui,2),ba(di,1),k)*0.5+0.5*f(ba(di,2),ba(ui,1),k);%f(ba(ui,2),ba(di,1),k);
end
% for k=kxv
%     xppl(k+kx+1)=f(ba(ui,1),ba(di,1),k);
%     xpnl(k+kx+1)=f(ba(ui,1),ba(di,2),k);
%     xnpl(k+kx+1)=f(ba(ui,2),ba(di,1),k);
% end
% for k=kxv
%     xppl(k+kx+1)=f(conj(ba(ui,1)),ba(ui,1),k)*0.5i-f(conj(ba(di,1)),ba(di,1),k)*0.5i+f(conj(ba(ui,1)),ba(di,1),k)*0.5+0.5*f(conj(ba(di,1)),ba(ui,1),k);%f(ba(ui,1),ba(di,1),k);
%     xpnl(k+kx+1)=f(conj(ba(ui,1)),ba(ui,2),k)*0.5i-f(conj(ba(di,1)),ba(di,2),k)*0.5i+f(conj(ba(ui,1)),ba(di,2),k)*0.5+0.5*f(conj(ba(di,1)),ba(ui,2),k);%f(ba(ui,1),ba(di,2),k);
%     xnpl(k+kx+1)=f(conj(ba(ui,2)),ba(ui,1),k)*0.5i-f(conj(ba(di,2)),ba(di,1),k)*0.5i+f(conj(ba(ui,2)),ba(di,1),k)*0.5+0.5*f(conj(ba(di,2)),ba(ui,1),k);%f(ba(ui,2),ba(di,1),k);
% end
Grel=kappa*sum(abs(xpnl).^2+abs(xnpl).^2);
Gdeph=0.5*kappa*sum(abs(xpnl).^2+abs(xnpl).^2+4*abs(xppl).^2);
rhoss=kappa*sum(abs(xnpl).^2)/Grel;
wxl=linspace(-30,30,1001)*kappa;
nx=length(wxl);
sw=zeros(1,nx);
% swa=zeros(1,nx);
for ii=1:nx
    sw(ii)=4*rhoss*(1-rhoss)*sum(abs(xppl).^2.*Grel./(Grel.^2+(wxl(ii)-kxv*delta).^2))...
        +sum(abs(xpnl).^2*rhoss*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta-Der).^2))...
        +sum(abs(xnpl).^2*(1-rhoss)*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta+Der).^2));
%     swa(ii)=sum((1-2*rhoss)*xpnl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2)+(2*rhoss-1)*xnpl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2));
end
figure(1)
swa=kappa*(0.25*0.5*kappa./(0.25*kappa^2+wxl.^2)+0.25*sum(0.75*kappa*repmat(besselj(kxv,2*A/delta).^2,length(wxl),1).'./(0.75^2*kappa^2+(repmat(wxl,length(kxv),1)-repmat(kxv*delta,length(wxl),1).').^2),1));
xx=wxl/kappa;
yy=delta*ones(1,nx);
% plot3(xx,yy,kappa*sw,'--'),hold on
plot(xx,sw),hold on
% save('sweq10kD0kr1delta0p01_offset_new.txt','xx','yy','swa','-ascii')