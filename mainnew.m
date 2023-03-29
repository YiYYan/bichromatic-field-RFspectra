% clc;
clear all;
global wx wz kx ky Nk Nh kxv kyv ui di Nkx Nky kappa kxy vd kxm weights vdm weighta
kx=200;% for the first mode
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
w0=1.0;wx=0.2;Ax=1;wz=1.0+wx;r=0.25;kappa=1.0;phi=pi/2;
Davg=0.0;
d1=wx/2+Davg;
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
vdm=repmat(vd,1,Nk);
uinv=inv(u);
tmpk=kappa./(1i*vd).*uinv(:,zeta0);
A1=uinv(:,1:3:end-2)*u(3:3:end,:)*tmpk;
A2=uinv(:,3:3:end)*u(2:3:end-1,:)*tmpk;
A3=(uinv(:,3:3:end).*(kappa./(1i*(vdm-kxm))))*u(2:3:end-1,:)*tmpk;
weights=ut(:,alpha0).*(0.5*uinv(:,alpha0)-0.5*A1+A2-A3);
weighta=ut(:,alpha0).*(A1-2*A2);
wxl=linspace(-40,40,1000)*kappa;
nx=length(wxl);
spfgft=zeros(1,nx);
spfgfta=zeros(1,nx);
for ii=1:nx
    spfgft(ii)=rfsfunnew(wxl(ii));
    spfgfta(ii)=rfafunnew(wxl(ii));
end
% figure(1);plot(wxl-wx/2,2*300*kappa*real(spfgft)),hold on
figure(2);plot(wxl-wx/2,kappa*real(spfgfta)),hold on