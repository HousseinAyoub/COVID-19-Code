function [Lambda]=ForceofInfection(x,beta,epsi,Nage,Ncomp)

eps=1e-10;

rhoN=zeros(Nage,1);

rhoN=sum(x(:,1:8),2);
rhoNtot=sum(rhoN(1:Nage));

%Mixing Matrix
%%%Age group
II=eye(Nage);
H=zeros(Nage,Nage);
for a=1:Nage
    H(:,a)=epsi.*II(:,a)+(1-epsi).*(rhoN(a)/(rhoNtot+eps)).*ones(Nage,1);
end

LambdaR=zeros(Nage,1);

for a=1:Nage 
    LambdaR(:)=LambdaR(:)+beta(:).*H(:,a).*((x(a,3)+x(a,4)+x(a,6))/(rhoN(a)+eps));
end
Lambda=LambdaR; 


