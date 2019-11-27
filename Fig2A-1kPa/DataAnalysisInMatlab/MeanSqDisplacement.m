addpath('/path/to/functions.')

Model='/path/to/model/';
variable='sqdis';
fname='Default';
folder = ['/path/to/InelasticityYoungsDiffusion/Data_Pars003-'];
pars=[7000 8000 9000 10000 11000 12000 14000 16000 18000 20000 25000 32000 ]/1000;
N=length(pars);
sim=100;
maxt=1000;
tstep=10;


sqdis=zeros(N,sim,length(0:tstep:maxt));
for y = 1:N
    y
    tic
    tindex=1;
    for tt=0:tstep:maxt
        tt
        for s=1:sim
        folder1 = [folder,num2str(y, '%03i'),'Sim',num2str(s, '%03i'),'/'];
        sqdis(y,s,tindex) = get_timeval(folder1,variable,tt);
         if(sqdis(y,s,tindex)==2^32)
             cd(Model)
             system(['./CPMFEM ',folder1])
             sqdis(y,s,tindex) = get_timeval(folder1,variable,tt);
         end
        end
        tindex=tindex+1;
    end
    toc
end

save(['/ufs/rens/CPM_FEM_Project/Results/ForPaperRene/DifferentLambdaValuesDiffusion/sqdis03.mat'],'sqdis');

%all saved




load('/ufs/rens/CPM_FEM_Project/Results/ForPaperRene/DifferentLambdaValuesDiffusion/sqdis02.mat')
meansqdis=squeeze(mean(sqdis,2));
vals=pars;


N=length(vals);
x=0:tstep:maxt;
M=meansqdis;
xl='time';
yl='ensemble MSD';
l='Youngs';

figure(1)
ColorSet=varycolor(6);
set(gca,'ColorOrder',ColorSet)
hold all;
for i=1:N 
    if(i<N/2+1)
plot(x,M(i,:),'LineWidth',2)
    else
plot(x,M(i,:),'--','LineWidth',2)
    end
end
xlabel(xl)
ylabel(yl)
h=legend('7000' ,'8000' ,'9000', '10000', '11000', '12000', '14000', '16000', '18000', '20000', '25000', '32000','Location','NorthWest')
v = get(h,'title');
set(v,'string','stiffness (Pa)');



tstart=600;
tstart=find(x==tstart);
tend=1000;
tend=find(x==tend);
dispcoeffs=[];
confints=[];
for y = 1:N
    ft = fittype('poly1');    
    [fo, gof] = fit((0:tstart)', meansqdis(y,1:tstart+1)',ft, 'Lower',[0 0],'Upper',[100 100]);
    ci = confint(fo);
    dispcoeffs(y,1)=fo.p1/4;
    confints(y,1:2)=ci(:,1)'/4;
    
    [fo, gof] = fit((tstart+1:tend)', meansqdis(y,tstart+1:tend)',ft,'Lower',[0 0],'Upper',[100 100]);
    ci = confint(fo);
    dispcoeffs(y,2)=fo.p1/4;
    confints(y,3:4)=ci(:,1)'/4;
end


figure(2)
mpar=dispcoeffs(:,1);
L=confints(:,1)-dispcoeffs(:,1);
U=confints(:,2)-dispcoeffs(:,1);
errorbar(1:N,mpar,L,U,'*black');
xlabel('stiffness (kPa)');
ylabel('D')
set(gca,'XTick',1:N)
xlim([0 N+1])
set(gca, 'XTickLabel', pars)

figure(3)
mpar=dispcoeffs(:,2);
L=confints(:,3)-dispcoeffs(:,2);
U=confints(:,4)-dispcoeffs(:,2);
errorbar(1:N,mpar,L,U,'*black');
xlabel('stiffness (kPa)');
ylabel('D')
set(gca,'XTick',1:N)
xlim([0 N+1])
set(gca, 'XTickLabel', pars)