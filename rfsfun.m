function rfs=rfsfun(x)
global vd ut uinv alpha0 A1 A2 A3
s1=ut(:,alpha0).*uinv(:,alpha0)*0.5./(1i*(x+vd));
s2=-0.5*ut(:,alpha0)./(1i*(x+vd)).*A1;
s2p=ut(:,alpha0)./(1i*(x+vd)).*A2;
s3=-ut(:,alpha0)./(1i*(x+vd)).*A3;
rfs=sum(s1+s2+s2p+s3);
end