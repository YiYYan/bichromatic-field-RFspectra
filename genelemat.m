function [Im,Dm,Cm]=genelemat(k)
Nk=2*k+1;
Im=eye(Nk);
Dm=diag(-k:k);
Cm=diag(ones(1,Nk-1),1);
end