%% COVID-19 Model
%Qatar University and WCMC-Q-Research Dept.
% Date: 18/07/2020
% Authors: Houssein H. Ayoub et al.
%%
clc
clear all
close all

load PopSize %China Demography including the CFR data
Demog='Demo.xlsx';  %%Demography of all countries as per the UN database
PopSizeAll= xlsread(Demog);

%% define the number of stages and INFECTION groups.
nst =10; % NUMBER OF STAGES (e.g. S, E, IM, IS, IC,...)
NAge =9; % number of age groups (10 years age band)

%% Time scale (days)
t0=0;         % Start time
tf=700;         % Stop time
dt=0.5;          % Time setp
tspan=t0:dt:tf;  % Timespan to use in the ode

eta=(ones(NAge,1)).*(1/(10*365)); %% Transition rate from one age group to the next age group
eta(NAge)=0;

%%%Fraction of Mild/Asym, Severe, and Critical cases, respectively
fM=[0.889 0.88 0.88 0.88 0.88 0.825 0.712 0.712 0.712];
fS=[0.111 0.099 0.099 0.099 0.099 0.103 0.078 0.078 0.078];
fC=[0 0.022 0.022 0.022 0.022 0.072 0.209 0.209 0.209];

%%%%The 500 model parameters were generated using Latin Hypercube sampling
load('deltan','deltan')
load('nuMn','nuMn')
load('nuSn','nuSn')
load('nuCn','nuCn')
load('nuSIDn','nuSIDn')
load('nuCIDn','nuCIDn')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Results of the 500 model fits (Please refer to: Characterizing key attributes of the epidemiology of COVID-19 in China: Model-based estimations)
load('EstimatedParamUncertainty500.mat','EstimatedParamUncertainty500')
zz=EstimatedParamUncertainty500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NCountry=159; %%Number of country included in the study

%%%Population size
N0=zeros(NAge,NCountry);
for Country=1:NCountry
mu=1/(PopSizeAll(Country,end)*365); %1/life expectancy
    for a=2:2:16
        N0(round(a/2),Country)=PopSizeAll(Country,a)+PopSizeAll(Country,a-1);  
    end
    N0(9,Country)=PopSizeAll(Country,17)+PopSizeAll(Country,18)+PopSizeAll(Country,19)+PopSizeAll(Country,20);
end
nruns=500;
tic
for Country=1:NCountry  %%%Parallel could be used here to minimize the time of runing
    for n=1:nruns
        Seed45=zz(20,n);
        monthsdist=round(zz(21,n));
        Seed45D=zz(22,n);
        monthsdistDeaths=round(zz(23,n));
        %%%%%%%%%%%%%%%%%%%%
        tinit=zz(1,n);
        a0=zz(2,n); %%contact rate
        a1=zz(3,n); %%contact rate
        b1=0; 
        c1=0;
        sigma=zeros(NAge,1);
        for j=1:NAge
            sigma(j)=zz(j+4,n); %Susceptibility by age group
        end
        FactorDeath=zz(14,n); 
        epsi=zz(15,n); %%Degree of assortativeness
        
        CFR=PopSize(1:NAge,6); %Case Fatality Rate by age group
        alpha=zeros(NAge,1); 
        for a=1:NAge
            alpha(a)=FactorDeath*CFR(a); %disease Mortality rate
        end
        
        %Initial conditions
        initprev=1/NAge; %initialnumber of COVID-19 cases
        x0L=zeros(NAge,nst);
        for a=1:NAge
        x0L(a,1)= N0(a,Country)-initprev;
        x0L(a,2)= initprev;
        end
        x0=zeros(nst*NAge,1);
        ii=1:nst;
        for j=1:NAge
            isj= ii+(j-1)*nst;
            x0(isj)=x0L(j,ii);
        end

        [T,x] = ode15s(@(t,x)COVID19_Model(t,x,a0,a1,b1,c1,sigma,alpha,deltan(n),NAge,mu,epsi,nst,eta,fM,fS,fC,nuMn(n),nuSn(n),nuCn(n),nuSIDn(n),nuCIDn(n)),tspan,x0);%,options);

        TT=length(T);

        SusA(:,:)=x(:,1:nst:nst*NAge); 
        LatentA(:,:)=x(:,2:nst:nst*NAge);

        InfectedIMA(:,:)=x(:,3:nst:nst*NAge);

        InfectedISA(:,:)=x(:,4:nst:nst*NAge);
        InfectedDSA(:,:)=x(:,5:nst:nst*NAge);

        InfectedICA(:,:)=x(:,6:nst:nst*NAge);
        InfectedDCA(:,:)=x(:,7:nst:nst*NAge);

        RecoveredA(:,:)=x(:,8:nst:nst*NAge);

        CumulativeIncidence(:,:)=x(:,9:nst:nst*NAge);

        CumulativeDeaths(:,:)=x(:,10:nst:nst*NAge);

        Total(:,:,n,Country)=SusA(:,:)+LatentA(:,:)+InfectedIMA(:,:)+(InfectedISA(:,:)+InfectedDSA(:,:))+(InfectedICA(:,:)+InfectedDCA(:,:))+RecoveredA(:,:);

        InfectedIS(:)=sum(x(:,4:nst:nst*NAge),2)';
        InfectedIC(:)=sum(x(:,6:nst:nst*NAge),2)';

        for t=1:TT
            for a=1:NAge
                IncidenceInfectedAge(t,a)=deltan(n)*LatentA(t,a);

                MortalityCasesTA(t,a)=alpha(a)*InfectedDCA(t,a);

                IncidenceMildA(t,a)=fM(a)*deltan(n)*LatentA(t,a);
            end
             IncidenceMild(t)=sum(IncidenceMildA(t,:));
             IncidenceMildT(t)=dt*sum(IncidenceMild(1:t));

            CumInfectedIS(t)=dt*nuSIDn(n)*sum(InfectedIS(1:t));
            CumInfectedIC(t)=dt*nuCIDn(n)*sum(InfectedIC(1:t));

            IncidenceInfected(t)=sum(IncidenceInfectedAge(t,:));
            MortalityCases(t)=sum(MortalityCasesTA(t,:));
            CumulativeIncidenceT(t)=sum(CumulativeIncidence(t,:));
            CumulativeDeathsT(t)=sum(CumulativeDeaths(t,:));
        end
        for t=t0:tf-1
            IncidenceInfectedT(t-t0+1,n,Country)=trapz(IncidenceInfected((t-t0)/dt+1:(t-t0+dt)/dt+1));
            MortalityCasesT(t-t0+1,n,Country)=trapz(MortalityCases((t-t0)/dt+1:(t-t0+dt)/dt+1));
            CumulativeIncidenceTT(t-t0+1,n,Country)=trapz(CumulativeIncidenceT((t-t0)/dt+1:(t-t0+dt)/dt+1));
            CumulativeDeathsTT(t-t0+1,n,Country)=trapz(CumulativeDeathsT((t-t0)/dt+1:(t-t0+dt)/dt+1));

            CumInfectedIST(t-t0+1,n,Country)=trapz(CumInfectedIS((t-t0)/dt+1:(t-t0+dt)/dt+1));
            CumInfectedICT(t-t0+1,n,Country)=trapz(CumInfectedIC((t-t0)/dt+1:(t-t0+dt)/dt+1));    
            IncidenceMildTT(t-t0+1,n,Country)=trapz(IncidenceMildT((t-t0)/dt+1:(t-t0+dt)/dt+1));

            for a=1:NAge
                Totalt(t-t0+1,a,n,Country)=trapz(Total((t-t0)/dt+1:(t-t0+dt)/dt+1,a,n,Country));
            end
        end
        
        for tx=1:length(tspan)
            beta(tx,n)=a0*(1+a1);
            for a=1:NAge
               R0a(tx,a,n,Country)=(Total(tx,a,n,Country)/sum(Total(tx,1:NAge,n,Country)))*fM(a)*((beta(tx,n)*deltan(n)*sigma(a))/((deltan(n)+mu+eta(a))*(nuMn(n)+mu+eta(a))))+(Total(tx,a,n,Country)/sum(Total(tx,1:NAge,n,Country)))*fS(a)*((beta(tx,n)*deltan(n)*sigma(a))/((deltan(n)+mu+eta(a))*(nuSIDn(n)+mu+eta(a))))+(Total(tx,a,n,Country)/sum(Total(tx,1:NAge,n,Country)))*fC(a)*((beta(tx,n)*deltan(n)*sigma(a))/((deltan(n)+mu+eta(a))*(nuCIDn(n)+mu+eta(a))));
            end
            R0(tx,n,Country)=sum(R0a(tx,1:NAge,n,Country));
        end
        [Va,PeakIncidence(n,Country)]=max(IncidenceInfectedT(:,n,Country));

        InfectiousCapita(n,Country)=100.*CumulativeIncidenceTT(end,n,Country)/sum(Totalt(end,:,n,Country),2);

        SevereCapita(n,Country)=100.*CumInfectedIST(end,n,Country)/sum(Totalt(end,:,n,Country),2);

        CriticalCapita(n,Country)=100.*CumInfectedICT(end,n,Country)/sum(Totalt(end,:,n,Country),2);

        DeathsCapita(n,Country)=100.*CumulativeDeathsTT(end,n,Country)/sum(Totalt(end,:,n,Country),2);

        MildCapita(n,Country)=100.*IncidenceMildTT(end,n,Country)/sum(Totalt(end,:,n,Country),2);
    end
end
toc

for Country=1:NCountry
   for n=1:nruns
       R0U{Country}(n)=R0(100,n,Country); 
       if (PeakIncidence(n,Country)>10) %%%Filter to remove countries with R0<1
           CumulativeIncidenceTTU{Country}(n)=CumulativeIncidenceTT(end,n,Country);
           CumulativeDeathsTU{Country}(n)=CumulativeDeathsTT(end,n,Country);
           
           PeakIncidenceU{Country}(n)=PeakIncidence(n,Country); 
           InfectiousCapitaU{Country}(n)=InfectiousCapita(n,Country); 
           SevereCapitaU{Country}(n)=SevereCapita(n,Country); 
           CriticalCapitaU{Country}(n)=CriticalCapita(n,Country);
           SevreCrticalCapitaU{Country}(n)=SevereCapita(n,Country)+CriticalCapita(n,Country);
           DeathsCapitaU{Country}(n)=DeathsCapita(n,Country); 
           MildCapitaU{Country}(n)=MildCapita(n,Country); 
       else
           PeakIncidenceU{Country}(n)=0; 
           InfectiousCapitaU{Country}(n)=0; 
           SevereCapitaU{Country}(n)=0; 
           CriticalCapitaU{Country}(n)=0;
           SevreCrticalCapitaU{Country}(n)=0;
           DeathsCapitaU{Country}(n)=0; 
           MildCapitaU{Country}(n)=0;
           
           CumulativeIncidenceTTU{Country}(n)=0;
           CumulativeDeathsTU{Country}(n)=0;
       end
   end
end
for Country=1:NCountry
    CumulativeIncidenceTTU{Country}(:,~any(CumulativeIncidenceTTU{Country},1))=[]; 
    CumulativeDeathsTU{Country}(:,~any(CumulativeDeathsTU{Country},1))=[]; 

    PeakIncidenceU{Country}(:,~any(PeakIncidenceU{Country},1))=[];  
    InfectiousCapitaU{Country}(:,~any(InfectiousCapitaU{Country},1))=[];  
    SevreCrticalCapitaU{Country}(:,~any(SevreCrticalCapitaU{Country},1))=[];  
    DeathsCapitaU{Country}(:,~any(DeathsCapitaU{Country},1))=[]; 
    MildCapitaU{Country}(:,~any(MildCapitaU{Country},1))=[]; 
end    

for Country=1:NCountry
    %%%%Most probable cumulative incidence
    [BinHeightCumulativeIncidenceTTU,BinCenterCumulativeIncidenceTTU]=createFit(CumulativeIncidenceTTU{Country}(:));
    [M,IICumulativeIncidenceTTU(Country)]=max(BinHeightCumulativeIncidenceTTU); 
    MostProbableCumulativeIncidenceTTU(Country)=BinCenterCumulativeIncidenceTTU(IICumulativeIncidenceTTU(Country));
end 
 for Country=1:NCountry
     %%%%Most probable cumulative deaths
    [BinHeightCumulativeDeathsTU,BinCenterCumulativeDeathsTU]=createFit(CumulativeDeathsTU{Country}(:));
    [M,IICumulativeDeathsTU(Country)]=max(BinHeightCumulativeDeathsTU); 
    MostProbableCumulativeDeathsTU(Country)=BinCenterCumulativeDeathsTU(IICumulativeDeathsTU(Country));
end

for Country=1:NCountry
    %%%%Most probable R0
    [BinHeightR0,BinCenterR0]=createFit(R0U{Country}(:));
    [M,IIR0(Country)]=max(BinHeightR0); 
    MostProbableR0(Country)=BinCenterR0(IIR0(Country)); 
    
    R0U1{Country}=sort(R0U{Country});
    R0U2{Country}=R0U1{Country}(round(length(R0U1{Country})*0.025):round(length(R0U1{Country})*(1-0.025)));
    LowBoundR0U(Country)=R0U2{Country}(1);
    UpperBoundR0U(Country)=R0U2{Country}(length(R0U2{Country}));
    MeanR0U(Country)=mean(R0U2{Country});
    MedianR0U(Country)=median(R0U2{Country});
end



for Country=1:NCountry    
     %%%%Most probable incidence peak
    [BinHeightPeakIncidenceU,BinCenterPeakIncidenceU]=createFit(PeakIncidenceU{Country}(:));
    [M,IIPeakIncidenceU(Country)]=max(BinHeightPeakIncidenceU); 
    MostProbablePeakIncidenceU(Country)=BinCenterPeakIncidenceU(IIPeakIncidenceU(Country)); 
    
    PeakIncidenceU1{Country}=sort(PeakIncidenceU{Country});
    PeakIncidenceU2{Country}=PeakIncidenceU1{Country}(round(length(PeakIncidenceU1{Country})*0.025):round(length(PeakIncidenceU1{Country})*(1-0.025)));
    LowBoundPeakIncidenceU(Country)=PeakIncidenceU2{Country}(1);
    UpperBoundPeakIncidenceU(Country)=PeakIncidenceU2{Country}(length(PeakIncidenceU2{Country}));
    MeanPeakIncidenceU(Country)=mean(PeakIncidenceU2{Country});
    MedianPeakIncidenceU(Country)=median(PeakIncidenceU2{Country});
end  



for Country=1:NCountry       
     %%%%Most probable Infections per capita
    [BinHeightInfectiousCapitaU,BinCenterInfectiousCapitaU]=createFit(InfectiousCapitaU{Country}(:));
    [M,IIInfectiousCapitaU(Country)]=max(BinHeightInfectiousCapitaU); 
    MostProbableInfectiousCapitaU(Country)=BinCenterInfectiousCapitaU(IIInfectiousCapitaU(Country)); 
    
    InfectiousCapitaU1{Country}=sort(InfectiousCapitaU{Country});
    InfectiousCapitaU2{Country}=InfectiousCapitaU1{Country}(round(length(InfectiousCapitaU1{Country})*0.025):round(length(InfectiousCapitaU1{Country})*(1-0.025)));
    LowBoundInfectiousCapitaU(Country)=InfectiousCapitaU2{Country}(1);
    UpperBoundInfectiousCapitaU(Country)=InfectiousCapitaU2{Country}(length(InfectiousCapitaU2{Country}));
    MeanInfectiousCapitaU(Country)=mean(InfectiousCapitaU2{Country});
    MedianInfectiousCapitaU(Country)=median(InfectiousCapitaU2{Country});
end 

for Country=1:NCountry
    %%%%Most probable Deaths per capita 
    [BinHeightDeathsCapitaU,BinCenterDeathsCapitaU]=createFit(DeathsCapitaU{Country}(:));
    [M,IIDeathsCapitaU(Country)]=max(BinHeightDeathsCapitaU); 
    MostProbableDeathsCapitaU(Country)=BinCenterDeathsCapitaU(IIDeathsCapitaU(Country)); 
    
    DeathsCapitaU1{Country}=sort(DeathsCapitaU{Country});
    DeathsCapitaU2{Country}=DeathsCapitaU1{Country}(round(length(DeathsCapitaU1{Country})*0.025):round(length(DeathsCapitaU1{Country})*(1-0.025)));
    LowBoundDeathsCapitaU(Country)=DeathsCapitaU2{Country}(1);
    UpperBoundDeathsCapitaU(Country)=DeathsCapitaU2{Country}(length(DeathsCapitaU2{Country}));
    MeanDeathsCapitaU(Country)=mean(DeathsCapitaU2{Country});
    MedianDeathsCapitaU(Country)=median(DeathsCapitaU2{Country});
end 

for Country=1:NCountry         
    %%%%Most probable mild per capita
    [BinHeightMildCapitaU,BinCenterMildCapitaU]=createFit(MildCapitaU{Country}(:));
    [M,IIMildCapitaU(Country)]=max(BinHeightMildCapitaU); 
    MostProbableMildCapitaU(Country)=BinCenterMildCapitaU(IIMildCapitaU(Country)); 
    
    MildCapitaU1{Country}=sort(MildCapitaU{Country});
    MildCapitaU2{Country}=MildCapitaU1{Country}(round(length(MildCapitaU1{Country})*0.025):round(length(MildCapitaU1{Country})*(1-0.025)));
    LowBoundMildCapitaU(Country)=MildCapitaU2{Country}(1);
    UpperBoundMildCapitaU(Country)=MildCapitaU2{Country}(length(MildCapitaU2{Country}));
    MeanMildCapitaU(Country)=mean(MildCapitaU2{Country});
    MedianMildCapitaU(Country)=median(MildCapitaU2{Country});
end

for Country=1:NCountry
      %%%%Most probable severe and critical per capita
    [BinHeightSevreCrticalCapitaU,BinCenterSevreCrticalCapitaU]=createFit(SevreCrticalCapitaU{Country}(:));
    [M,IISevreCrticalCapitaU(Country)]=max(BinHeightSevreCrticalCapitaU); 
    MostProbableSevreCrticalCapitaU(Country)=BinCenterSevreCrticalCapitaU(IISevreCrticalCapitaU(Country)); 
    
    SevereCapitaU1{Country}=sort(SevreCrticalCapitaU{Country});
    SevereCapitaU2{Country}=SevereCapitaU1{Country}(round(length(SevereCapitaU1{Country})*0.025):round(length(SevereCapitaU1{Country})*(1-0.025)));
    LowBoundSevereCapitaU(Country)=SevereCapitaU2{Country}(1);
    UpperBoundSevereCapitaU(Country)=SevereCapitaU2{Country}(length(SevereCapitaU2{Country}));
    MeanSevereCapitaU(Country)=mean(SevereCapitaU2{Country});
    MedianSevereCapitaU(Country)=median(SevereCapitaU2{Country});
end

%%%%%%%%%%%%%%Figure 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%AFRO
xlim([1,160])
bar(1:43,MostProbableR0(1:43),'b')
hold on
%%%EMRO
bar(44:64,MostProbableR0(69:89),'g')
%%%SEARO
bar(65:73,MostProbableR0(136:144),'b')
%%%AMRO
bar(74:98,MostProbableR0(44:68),'r')
%%%WPRO
bar(99:113,MostProbableR0(145:159),'k')
%%%EURO
bar(114:159,MostProbableR0(90:135),'y')

hold on
plot([0.5,160],[1,1],'r','LineWidth',3)

median(MostProbableR0(1:43))
median(MostProbableR0(69:89))
median(MostProbableR0(136:144))
median(MostProbableR0(44:68))
median(MostProbableR0(145:159))
median(MostProbableR0(90:135))



SelectedCountryMostProbableR0=[MostProbableR0(111) MostProbableR0(66) MostProbableR0(147) MostProbableR0(46) MostProbableR0(139) MostProbableR0(138) MostProbableR0(71) MostProbableR0(80) MostProbableR0(32)];
SelectedCountryMostProbablePeakIncidence=[MostProbablePeakIncidenceU(111) MostProbablePeakIncidenceU(66) MostProbablePeakIncidenceU(147) MostProbablePeakIncidenceU(46) MostProbablePeakIncidenceU(139) MostProbablePeakIncidenceU(138) MostProbablePeakIncidenceU(71) MostProbablePeakIncidenceU(80) MostProbablePeakIncidenceU(32)];
SelectedCountryMostProbableInfectiousCapitaU=[MostProbableInfectiousCapitaU(111) MostProbableInfectiousCapitaU(66) MostProbableInfectiousCapitaU(147) MostProbableInfectiousCapitaU(46) MostProbableInfectiousCapitaU(139) MostProbableInfectiousCapitaU(138) MostProbableInfectiousCapitaU(71) MostProbableInfectiousCapitaU(80) MostProbableInfectiousCapitaU(32)];
SelectedCountryMostProbableDeathsCapitaU=[MostProbableDeathsCapitaU(111) MostProbableDeathsCapitaU(66) MostProbableDeathsCapitaU(147) MostProbableDeathsCapitaU(46) MostProbableDeathsCapitaU(139) MostProbableDeathsCapitaU(138) MostProbableDeathsCapitaU(71) MostProbableDeathsCapitaU(80) MostProbableDeathsCapitaU(32)];
SelectedCountryMostProbableMildCapitaU=[MostProbableMildCapitaU(111) MostProbableMildCapitaU(66) MostProbableMildCapitaU(147) MostProbableMildCapitaU(46) MostProbableMildCapitaU(139) MostProbableMildCapitaU(138) MostProbableMildCapitaU(71) MostProbableMildCapitaU(80) MostProbableMildCapitaU(32)];
SelectedCountryMostProbableSevreCrticalCapitaU=[MostProbableSevreCrticalCapitaU(111) MostProbableSevreCrticalCapitaU(66) MostProbableSevreCrticalCapitaU(147) MostProbableSevreCrticalCapitaU(46) MostProbableSevreCrticalCapitaU(139) MostProbableSevreCrticalCapitaU(138) MostProbableSevreCrticalCapitaU(71) MostProbableSevreCrticalCapitaU(80) MostProbableSevreCrticalCapitaU(32)];

%%%%%%%%%%%%%%Figure 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(211)
bar(SelectedCountryMostProbableR0);
set(gca,'XTick',[1 2 3 4 5 6 7 8 9],'XTickLabel',{'Italy','US','China','Brazil','Indonesia','India','Egypt','Pakistan','Niger'}); 
xlim([0.5 9.5])
text(1:length(SelectedCountryMostProbableR0),round(SelectedCountryMostProbableR0,2),num2str(round(SelectedCountryMostProbableR0,2)'),'vert','bottom','horiz','center','FontWeight','bold'); 
hold on
plot([0.5,9.5],[1,1],'r','LineWidth',2)
ylabel('Basic reproduction number (R_{0})')
LowerL=SelectedCountryMostProbableR0-[2.02	1.76	1.75	1.57	1.41	1.37	1.22	1.14	0.83];
UpperL=[3.27	2.86	2.86	2.57	2.32	2.25	2.01	1.88	1.37]-SelectedCountryMostProbableR0; 
hh2011=errorbar(1:length(SelectedCountryMostProbableR0),SelectedCountryMostProbableR0,LowerL,UpperL,'MarkerSize',5,'Marker','*',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);
box off
% 
subplot(212)
%Peak
hold on
bar(SelectedCountryMostProbablePeakIncidence(1:9-1));
set(gca,'XTick',[1 2 3 4 5 6 7 8 9],'XTickLabel',{'Italy','US','China','Brazil','Indonesia','India','Egypt','Pakistan','Niger'}); 
xlim([0.5 9.5])
text(1:length(SelectedCountryMostProbablePeakIncidence),round(SelectedCountryMostProbablePeakIncidence),num2str(round(SelectedCountryMostProbablePeakIncidence)'),'vert','bottom','horiz','center','FontWeight','bold'); 
ylabel('Day at peak incidence')
LowerL=SelectedCountryMostProbablePeakIncidence(1:8)-[116	156	169	180	211	240	255	303];
UpperL=[124	160	172	190	242	283	350	479]-SelectedCountryMostProbablePeakIncidence(1:8); 
hh2011=errorbar(1:length(SelectedCountryMostProbableR0)-1,SelectedCountryMostProbablePeakIncidence(1:8),LowerL,UpperL,'MarkerSize',5,'Marker','*',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);
box off

%%%%%%%%%%%%%%Figure 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(211)
hold on
bar(SelectedCountryMostProbableInfectiousCapitaU(1:9-1));
set(gca,'XTick',[1 2 3 4 5 6 7 8 9],'XTickLabel',{'Italy','US','China','Brazil','Indonesia','India','Egypt','Pakistan','Niger'}); 
xlim([0.5 9.5])
text(1:length(SelectedCountryMostProbableInfectiousCapitaU),round(SelectedCountryMostProbableInfectiousCapitaU,1),num2str(round(SelectedCountryMostProbableInfectiousCapitaU,1)'),'vert','bottom','horiz','center','FontWeight','bold'); 
ylabel('Infections per capita (per 100 persons)')
LowerL=SelectedCountryMostProbableInfectiousCapitaU(1:8)-[69.6	57.9	58.7	49.3	39.9	37.7	26.7	21.8];
UpperL=[79.5	70.2	71.4	64.3	56.6	55.0	45.5	41.3]-SelectedCountryMostProbableInfectiousCapitaU(1:8); 
hh2011=errorbar(1:length(SelectedCountryMostProbableR0)-1,SelectedCountryMostProbableInfectiousCapitaU(1:8),LowerL,UpperL,'MarkerSize',5,'Marker','*',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);
box off
subplot(212)
hold on
bar(SelectedCountryMostProbableDeathsCapitaU(1:9-1));
set(gca,'XTick',[1 2 3 4 5 6 7 8 9],'XTickLabel',{'Italy','US','China','Brazil','Indonesia','India','Egypt','Pakistan','Niger'}); 
xlim([0.5 9.5])
text(1:length(SelectedCountryMostProbableDeathsCapitaU),round(SelectedCountryMostProbableDeathsCapitaU,1),num2str(round(SelectedCountryMostProbableDeathsCapitaU,1)'),'vert','bottom','horiz','center','FontWeight','bold'); 
ylabel('Deaths per capita (per 100 persons)')

LowerL=SelectedCountryMostProbableDeathsCapitaU(1:8)-[4.5	3.2	2.5	2.0	1.3	1.3	0.9	0.7];
UpperL=[5.1	3.8	3.0	2.4	1.7	1.7	1.3	1.1]-SelectedCountryMostProbableDeathsCapitaU(1:8); 
hh2011=errorbar(1:length(SelectedCountryMostProbableR0)-1,SelectedCountryMostProbableDeathsCapitaU(1:8),LowerL,UpperL,'MarkerSize',5,'Marker','*',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);
box off

%%%%%%%%%%%%%%Figure 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%MildCapita
figure
subplot(211)
hold on
bar(SelectedCountryMostProbableMildCapitaU(1:9-1));
set(gca,'XTick',[1 2 3 4 5 6 7 8 9],'XTickLabel',{'Italy','US','China','Brazil','Indonesia','India','Egypt','Pakistan','Niger'}); 
xlim([0.5 9.5])
text(1:length(SelectedCountryMostProbableMildCapitaU),round(SelectedCountryMostProbableMildCapitaU,1),num2str(round(SelectedCountryMostProbableMildCapitaU,1)'),'vert','bottom','horiz','center','FontWeight','bold'); 
ylabel('COVID-19 mild infections per capita (per 100 persons)')
LowerL=SelectedCountryMostProbableMildCapitaU(1:8)-[55.9	47.1	48.4	41.0	33.4	31.6	22.4	18.3];
UpperL=[64.1	57.4	59.1	53.7	47.6	46.3	38.4	35.0]-SelectedCountryMostProbableMildCapitaU(1:8); 
hh2011=errorbar(1:length(SelectedCountryMostProbableR0)-1,SelectedCountryMostProbableMildCapitaU(1:8),LowerL,UpperL,'MarkerSize',5,'Marker','*',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);
box off

subplot(212)
hold on
bar(SelectedCountryMostProbableSevreCrticalCapitaU(1:9-1));
set(gca,'XTick',[1 2 3 4 5 6 7 8 9],'XTickLabel',{'Italy','US','China','Brazil','Indonesia','India','Egypt','Pakistan','Niger'}); 
xlim([0.5 9.5])
text(1:length(SelectedCountryMostProbableSevreCrticalCapitaU),round(SelectedCountryMostProbableSevreCrticalCapitaU,1),num2str(round(SelectedCountryMostProbableSevreCrticalCapitaU,1)'),'vert','bottom','horiz','center','FontWeight','bold'); 
ylabel('COVID-19 severe and critical disease per capita (per 100 persons)')
LowerL=SelectedCountryMostProbableSevreCrticalCapitaU(1:8)-[13.7	10.9	10.3	8.4	6.5	6.2	4.3	3.4];
UpperL=[15.4	12.9	12.3	10.6	9.0	8.7	7.1	6.3]-SelectedCountryMostProbableSevreCrticalCapitaU(1:8); 
hh2011=errorbar(1:length(SelectedCountryMostProbableR0)-1,SelectedCountryMostProbableSevreCrticalCapitaU(1:8),LowerL,UpperL,'MarkerSize',5,'Marker','*',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[1 0 0]);
box off
