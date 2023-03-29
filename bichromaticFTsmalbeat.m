% clc;
clear all;
global wx wz kx ky Nk Nh kxv kyv ui di Nkx Nky kappa kxy vd kxm uinvzeta0 uinvalpha0 vdm ut z0 uinv A1 A2 A3 alpha0 zeta0
kx=400;% for the first mode
ky=0;% for the second mode, if consider single mode drive, set ky=0; 
kxv=-kx:kx;kyv=-ky:ky;
Nkx=length(kxv);Nky=length(kyv);
Nk=Nkx*Nky;
Nh=3*Nkx*Nky;
ui=1:3:Nh-2;
mi=2:3:Nh-1;
di=3:3:Nh;
[Imx,Dmx,Cmx]=genelemat(kx);
[Imy,Dmy,Cmy]=genelemat(ky);
%%=================在此处修改参数========================
%%wx=(w2-w1);d1=w0-w1;Ax=2*Omega
w0=1.0;wx=0.02;Ax=20;wz=1.0+wx;r=1;kappa=1.0;phi=pi/2;
Davg=0;
d1=Davg+wx/2;
% ws=0.5*(wx+wz);
%wx==w1;wz==w2;
%=========================================
Dd=[-d1-1i*0.5*kappa 0 Ax/4
   0 d1-1i*0.5*kappa -Ax/4
   Ax/2 -Ax/2 -1i*kappa];
Cxp=[0 0 0.25*r*Ax
     0 0 0
     0 -0.5*r*Ax 0];
Cxn=r*[0 0 0
       0 0 -0.25*Ax
       0.5*Ax 0 0];
% Wz=kron(Imx,wz*eye(3));
% Cc=kron(Imx,Cx2);
alpha0=find(kron(kyv==0,kron(kxv==0,[1,0,0]))==1);
beta0=find(kron(kyv==0,kron(kxv==0,[0,1,0]))==1);
zeta0=find(kron(kyv==0,kron(kxv==0,[0,0,1]))==1);
kxy=kron(ones(1,Nky),kxv)*wx+kron(kyv,ones(1,Nkx))*wz;
kxm=repmat(kxy,Nh,1);
Wx=wx*eye(3);
fmx=kron(Imx,Dd)+kron(Dmx,Wx)+kron(Cmx,Cxn)+kron(Cmx.',Cxp);
% fm=kron(Imy,fmx)+kron(Dmy,Wz)+kron(Cmy,Cc)+kron(Cmy.',conj(Cc));
[u,v]=eig(fmx);
ut=u.';
vd=diag(v);
uinv=inv(u);
vdm=repmat(vd,1,Nk);
uinvzeta0=repmat(uinv(:,zeta0),1,Nk);
uinvalpha0=repmat(uinv(:,alpha0),1,Nk);
ffm=repmat(uinv(:,zeta0)*kappa./(1i*vd),1,Nk);
A1=sum(uinv(:,1:3:Nh-2).*repmat(sum(ut(:,3:3:Nh).*ffm),Nh,1),2);
A2=sum(uinv(:,3:3:Nh).*repmat(sum(ut(:,2:3:Nh-1).*ffm),Nh,1),2);
A3=sum(uinv(:,3:3:Nh)*kappa./(1i*(vdm-kxm)).*repmat(sum(ut(:,2:3:Nh-1).*ffm),Nh,1),2);
z0=-1;
wxl=linspace(-30,30,1001)*kappa+wx/2;
nx=length(wxl);
spfgft=zeros(1,nx);
for ii=1:nx
    spfgft(ii)=rfsfun(wxl(ii));
end
swfl=kappa*real(spfgft);
xx=wxl-wx/2;
% yy=wx/2*ones(1,nx);
plot3(xx,yy,swfl),hold on
% plot(xx,swfl),hold on
save('swfl10kD0kr1delta0p01_offset.txt','xx','yy','swfl','-ascii')