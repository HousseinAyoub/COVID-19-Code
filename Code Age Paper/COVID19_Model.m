function [dxL,lambda]=COVID19_Model(t,x,a0,a1,b1,c1,sigma,alpha,delta,NAge,mu,epsi,nst,eta,fM,fS,fC,nuM,nuS,nuC,nuSID,nuCID)

xt=zeros(NAge,nst);
ii=1:nst;
for j=1:NAge
    isj= ii+(j-1)*nst;
    xt(j,ii)=x(isj);
end
beta=zeros(NAge,1);

for j=1:NAge
    %beta(j)=a0*(1+(a1/(1+exp((t-b1)/c1)))).*sigma(j); Wood Saxon was
    %assumed for China to account for the social distancing measures
    
    beta(j)=a0*(1+a1).*sigma(j); %%Worst-case scenario: There is no social distancing measures
end


lambda=ForceofInfection(xt,beta,epsi,NAge,8);
%%
%
%Age 0-9
dx=zeros(NAge,nst);
dx(1,1)= sum((alpha.*xt(:,7)),1)  - (mu+lambda(1)+eta(1))*xt(1,1);             %(S)
dx(1,2)= lambda(1)*xt(1,1) - (mu+delta+eta(1))*xt(1,2) ;       %(E)

dx(1,3)=fM(1)*delta*xt(1,2) - (mu+nuM+eta(1))*xt(1,3);      %(IM)

dx(1,4)=fS(1)*delta*xt(1,2) - (mu+nuSID+eta(1))*xt(1,4);      %(IS) 
dx(1,5)=nuSID*xt(1,4) - (mu+nuS+eta(1))*xt(1,5);      %(DS)

dx(1,6)=fC(1)*delta*xt(1,2) - (mu+nuCID+eta(1))*xt(1,6);      %(IC)
dx(1,7)=nuCID*xt(1,6) - (mu+nuC+eta(1)+alpha(1))*xt(1,7);      %(DC)

dx(1,8)= nuM*xt(1,3)+nuS*xt(1,5)+nuC*xt(1,7)-(mu+eta(1))*xt(1,8);          %(R)

dx(1,9)=delta*xt(1,2);    %%Cumulative incidence 
dx(1,10)=alpha(1)*xt(1,7);    %%Cumulative deaths  

%Age >9
for rst=2:NAge
dx(rst,1)=eta(rst-1)*xt(rst-1,1)  - (mu+lambda(rst)+eta(rst))*xt(rst,1);             %(S)
dx(rst,2)= eta(rst-1)*xt(rst-1,2)+ lambda(rst)*xt(rst,1) - (mu+delta+eta(rst))*xt(rst,2) ;       %(E)

dx(rst,3)=eta(rst-1)*xt(rst-1,3)+ fM(rst)*delta*xt(rst,2) - (mu+nuM+eta(rst))*xt(rst,3);      %(IM)

dx(rst,4)=eta(rst-1)*xt(rst-1,4)+ fS(rst)*delta*xt(rst,2) - (mu+nuSID+eta(rst))*xt(rst,4);      %(IS) 
dx(rst,5)=eta(rst-1)*xt(rst-1,5)+ nuSID*xt(rst,4) - (mu+nuS+eta(rst))*xt(rst,5);      %(DS)

dx(rst,6)=eta(rst-1)*xt(rst-1,6)+ fC(rst)*delta*xt(rst,2) - (mu+nuCID+eta(rst))*xt(rst,6);      %(IC)
dx(rst,7)=eta(rst-1)*xt(rst-1,7)+ nuCID*xt(rst,6) - (mu+nuC+eta(rst)+alpha(rst))*xt(rst,7);      %(DC)

dx(rst,8)=eta(rst-1)*xt(rst-1,8)+ nuM*xt(rst,3)+nuS*xt(rst,5)+nuC*xt(rst,7)-(mu+eta(rst))*xt(rst,8);          %(R)

dx(rst,9)=delta*xt(rst,2);    %%Cumulative incidence
dx(rst,10)=alpha(rst)*xt(rst,7);    %%Cumulative deaths  

end


dxL=zeros(NAge*nst,1);
ii=1:nst;
for j=1:NAge
    isj=ii+(j-1)*nst;
dxL(isj)= dx(j,ii);
end

end



