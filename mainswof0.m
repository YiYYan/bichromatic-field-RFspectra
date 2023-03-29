clc;clear;
A=10;kappa=1.0;
kxv=-25000:25000;
xppl=0.5*(kxv==0);
xnnl=-xppl;
wxl=0.0*A;%linspace(-40,40,1001)*kappa;
dxl=linspace(0.001,10.0,500);
nx=length(dxl);
sw=zeros(1,nx);
swa=zeros(1,nx);
for ii=1:nx
    delta=dxl(ii);
    xpnl=-0.5*besselj(kxv,2*A/delta);
    xnpl=0.5*besselj(-kxv,2*A/delta);
    Grel=kappa*sum(xpnl.^2+xnpl.^2);
    Gdeph=0.5*kappa*sum(xpnl.^2+xnpl.^2+4*xppl.^2);
    rhoss=kappa*sum(xnpl.^2)/Grel;
    sw(ii)=4*rhoss*(1-rhoss)*sum(xppl.^2.*Grel./(Grel.^2+(wxl-kxv*delta).^2))...
        +sum(xpnl.^2*rhoss*Gdeph./(Gdeph.^2+(wxl-kxv*delta).^2))...
        +sum(xnpl.^2*(1-rhoss)*Gdeph./(Gdeph.^2+(wxl-kxv*delta).^2));
    swa(ii)=kappa*(0.25*0.5*kappa./(0.25*kappa^2+wxl.^2)+0.25*sum(0.75*kappa*repmat(besselj(kxv,2*A/delta).^2,length(wxl),1).'./(0.75^2*kappa^2+(repmat(wxl,length(kxv),1)-repmat(kxv*delta,length(wxl),1).').^2),1));%sum((1-2*rhoss)*xpnl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2)+(2*rhoss-1)*xnpl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2));
end
figure(1)
% swa=kappa*(0.25*0.5*kappa./(0.25*kappa^2+wxl.^2)+0.25*sum(0.75*kappa*repmat(besselj(kxv,2*A/delta).^2,length(wxl),1).'./(0.75^2*kappa^2+(repmat(wxl,length(kxv),1)-repmat(kxv*delta,length(wxl),1).').^2),1));
xx=dxl;
plot(xx,sw,xx,swa),hold on
save('sweq10kD0kr1_sw0.txt','xx','swa','-ascii')