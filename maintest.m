A=0.01;kappa=0.001;delta=1e-6;
kxv=-25000:25000;
xppl=0.5*(kxv==0);
xnnl=-xppl;
xpnl=-0.5*besselj(kxv,2*A/delta);
xnpl=0.5*besselj(-kxv,2*A/delta);
Grel=kappa*sum(xpnl.^2+xnpl.^2);
Gdeph=0.5*kappa*sum(xpnl.^2+xnpl.^2+4*xppl.^2);
rhoss=kappa*sum(xnpl.^2)/Grel;
wxl=linspace(-40,40,1001)*kappa;
nx=length(wxl);
sw=zeros(1,nx);
swa=zeros(1,nx);
for ii=1:nx
    sw(ii)=4*rhoss*(1-rhoss)*sum(xppl.^2.*Grel./(Grel.^2+(wxl(ii)-kxv*delta).^2))...
        +sum(xpnl.^2*rhoss*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2))...
        +sum(xnpl.^2*(1-rhoss)*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2));
%     swa(ii)=sum((1-2*rhoss)*xpnl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2)+(2*rhoss-1)*xnpl.^2.*Gdeph./(Gdeph.^2+(wxl(ii)-kxv*delta).^2));
end
figure(1)
swa=kappa*(0.25*0.5*kappa./(0.25*kappa^2+wxl.^2)+0.25*sum(0.75*kappa*repmat(besselj(kxv,2*A/delta).^2,length(wxl),1).'./(0.75^2*kappa^2+(repmat(wxl,length(kxv),1)-repmat(kxv*delta,length(wxl),1).').^2),1));
xx=wxl/kappa;
plot(xx,kappa*sw,xx,swa),hold on
% save('sweq15k.txt','xx','swa','-ascii')
% figure(2)
% plot(wxl,kappa*swa,'-'),hold on 0.25*0.5*kappa./(0.25*kappa^2+wxl.^2)+