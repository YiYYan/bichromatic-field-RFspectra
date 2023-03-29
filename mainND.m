clc;clear;
sx=[1 0;0 -1];
sz=[0 1;1 0];
sy=-[0 -1i;1i 0];
sp=(sx+1i*sy)/2;
sn=(sx-1i*sy)/2;
kappa=1.0;
Davg=0;
delta=0.01*kappa;
A=10.0*kappa;r=1;
nn=1000;
KK=2*nn;
kxv=-KK:KK;
Nk=length(kxv);
Nh=2*Nk;
ui=1:2:Nh-1;
di=2:2:Nh;
z=A*(1+r)/delta;
jn1=@(n) besselj(2*n,z)*Davg-0.5*A*(1-r)*(besselj(2*n-1,z)-besselj(2*n+1,z));
jn2=@(n) besselj(2*n+1,z)*Davg-0.5*A*(1-r)*(besselj(2*n,z)-besselj(2*n+2,z));
H0=0.5*(Davg*besselj(0,z)+A*(1-r)*besselj(1,z))*sz;
isyp=[0 -1;1 0];
isyn=[0 1;-1 0];
tmp=zeros(Nh,Nh);
for n=1:nn
    cn1=diag(ones(1,Nk-2*n),2*n);
    cn2=diag(ones(1,Nk-2*(n-1)-1),2*(n-1)+1);
    tmp1=kron(cn1,0.5*jn1(n)*sz);
    tmp2=kron(cn2,0.5*jn2(n-1)*isyp);
    tmp=tmp+tmp1+tmp2;
end
HF=kron(eye(Nk),H0)+kron(diag(kxv),delta*eye(2))+tmp+tmp';
[u,v]=eig(HF);
rc=find(diag(v)<0.5*delta&diag(v)>=-0.5*delta);
er=diag(v(rc,rc));
ba=u(:,rc);
Der=er(1)-er(2);
f=@(a,b,k) sum(a(1-k*(k<0):Nk-k*(k>0)).*b(1+k*(k>0):Nk+k*(k<0)));
xppl=zeros(1,Nk);
xpnl=zeros(1,Nk);
xnpl=zeros(1,Nk);
% for k=kxv
%     xppl(k+kx+1)=f(ba(ui,1),ba(ui,1),k)*0.5-f(ba(di,1),ba(di,1),k)*0.5-f(ba(ui,1),ba(di,1),k)*0.5+0.5*f(ba(di,1),ba(ui,1),k);%f(ba(ui,1),ba(di,1),k);
%     xpnl(k+kx+1)=f(ba(ui,1),ba(ui,2),k)*0.5-f(ba(di,1),ba(di,2),k)*0.5-f(ba(ui,1),ba(di,2),k)*0.5+0.5*f(ba(di,1),ba(ui,2),k);%f(ba(ui,1),ba(di,2),k);
%     xnpl(k+kx+1)=f(ba(ui,2),ba(ui,1),k)*0.5-f(ba(di,2),ba(di,1),k)*0.5-f(ba(ui,2),ba(di,1),k)*0.5+0.5*f(ba(di,2),ba(ui,1),k);%f(ba(ui,2),ba(di,1),k);
% end
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
bau1=ba(ui,1);bau2=ba(ui,2);bad1=ba(di,1);bad2=ba(di,2);
bsel=besselj(kxv.',z);
bselr=bsel(end:-1:1);
% for k=kxv
%     xppl(k+KK+1)=0.5*f(bau1,bau1,k)-0.5*f(bad1,bad1,k)-0.5*(sum(bsel.*bau1)*bad1(k+KK+1)+sum(bau1.*bad1)*besselj(k,z)-besselj(k,z)*bau1(k+KK+1)*bad1(k+KK+1))...
%         +0.5*(sum(bselr.*bad1)*bau1(k+KK+1)+sum(bau1.*bad1)*besselj(-k,z)-besselj(-k,z)*bad1(k+KK+1)*bau1(k+KK+1));
%     xpnl(k+KK+1)=0.5*f(bau1,bau2,k)-0.5*f(bad1,bad2,k)-0.5*(sum(bsel.*bau1)*bad2(k+KK+1)+sum(bau1.*bad2)*besselj(k,z)-besselj(k,z)*bau1(k+KK+1)*bad2(k+KK+1))...
%         +0.5*(sum(bselr.*bad1)*bau2(k+KK+1)+sum(bad1.*bau2)*besselj(-k,z)-besselj(-k,z)*bad1(k+KK+1)*bau2(k+KK+1));
%     xnpl(k+KK+1)=0.5*f(bau2,bau1,k)-0.5*f(bad2,bad1,k)-0.5*(sum(bsel.*bau2)*bad1(k+KK+1)+sum(bau2.*bad1)*besselj(k,z)-besselj(k,z)*bau2(k+KK+1)*bad1(k+KK+1))...
%         +0.5*(sum(bselr.*bad2)*bau1(k+KK+1)+sum(bad2.*bau1)*besselj(-k,z)-besselj(-k,z)*bad2(k+KK+1)*bau1(k+KK+1));
% end
ba11=zeros(Nk,1);
ba12=zeros(Nk,1);
ba21=zeros(Nk,1);
ab11=zeros(Nk,1);
ab12=zeros(Nk,1);
ab21=zeros(Nk,1);
for k=kxv
    ba11(k+KK+1)=f(bad1,bau1,k);
    ba12(k+KK+1)=f(bad1,bau2,k);
    ba21(k+KK+1)=f(bad2,bau1,k);
    ab11(k+KK+1)=f(bau1,bad1,k);
    ab12(k+KK+1)=f(bau1,bad2,k);
    ab21(k+KK+1)=f(bau2,bad1,k);
end
for k=kxv
    xppl(k+KK+1)=0.5*f(bau1,bau1,k)-0.5*f(bad1,bad1,k)-0.5*f(bselr,ab11,k)+0.5*f(bsel,ba11,k);
    xpnl(k+KK+1)=0.5*f(bau1,bau2,k)-0.5*f(bad1,bad2,k)-0.5*f(bselr,ab12,k)+0.5*f(bsel,ba12,k);
    xnpl(k+KK+1)=0.5*f(bau2,bau1,k)-0.5*f(bad2,bad1,k)-0.5*f(bselr,ab21,k)+0.5*f(bsel,ba21,k);
end
Grel=kappa*sum(abs(xpnl).^2+abs(xnpl).^2);
Gdeph=0.5*kappa*sum(abs(xpnl).^2+abs(xnpl).^2+4*abs(xppl).^2);
rhoss=kappa*sum(abs(xnpl).^2)/Grel;
wxl=linspace(-50,50,1001)*kappa;
nx=length(wxl);
sw=zeros(1,nx);
% swa=zeros(1,nx);
for ii=1:nx
    sw(ii)=4*rhoss*(1-rhoss)*sum(abs(xppl).^2.*Grel./(Grel.^2+(wxl(ii)-kxv*delta).^2))...
        +sum(abs(xpnl).^2*rhoss*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta-Der).^2))...
        +sum(abs(xnpl).^2*(1-rhoss)*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta+Der).^2));
    %swa(ii)=sum((1-2*rhoss)*xpnl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2)+(2*rhoss-1)*xnpl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2));
end
figure(1)
% swa=kappa*(0.25*0.5*kappa./(0.25*kappa^2+wxl.^2)+0.25*sum(0.75*kappa*repmat(besselj(kxv,2*A/delta).^2,length(wxl),1).'./(0.75^2*kappa^2+(repmat(wxl,length(kxv),1)-repmat(kxv*delta,length(wxl),1).').^2),1));
xx=wxl/kappa;
yy=delta*ones(1,nx);
plot3(xx,yy,kappa*sw,'--'),hold on
% plot(xx,sw),hold on