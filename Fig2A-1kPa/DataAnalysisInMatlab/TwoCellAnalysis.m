addpath('/path/to/folder/')
folderstart = 'path/to/JYoungsTwoCells/Data_Pars';
Model='/path/to/model/';
savefolder = '/path/to/savefolder/';

parnrs=1:7;
pars=[500 1000 2000 4000 8000 10000 12000 14000 16000 32000]/1000;
N=length(pars);
sim=100;

checkexistfile='cpmfem00000.png'

MCS1=20;
MCS2=500;
tstep=1;
loopvec=MCS1:tstep:MCS2;
loopvecsqdis=0:1:5000;

for p = parnrs
    
parnr=parnrs(p);
folder = [folderstart,num2str(parnr, '%03i'),'-'];



contactcounts=zeros(N,sim);
contactlengths=zeros(N,sim,length(loopvec));
acute=zeros(N,sim);
sqdis1=zeros(N,sim,length(loopvecsqdis));
sqdis2=zeros(N,sim,length(loopvecsqdis));

for n = 1:N
    n
    tic
    for s = 1:sim
        folder1 = [folder,num2str(n, '%03i'),'Sim',num2str(s, '%03i'),'/'];
        if(~exist([folder1,checkexistfile]))
            cd(Model)
            system(['./CPMFEM ',folder1])
        end
        [c1,c2]=GetContacts(folder1,loopvec); 
        contactcounts(n,s)=c1;
        contactlengths(n,s,1:c1)=c2;
        acute(n,s)=GetObtuseAcute(folder1,loopvec);
        sqdis=loadvariableMCS([folder1,'sqdis.txt'],loopvecsqdis,2);
        sqdis1(n,s,:)=sqdis(:,1);
        sqdis2(n,s,:)=sqdis(:,2);

    end
    
    toc
end

save([savefolder,'acute-contact-sqdis',num2str(parnr, '%02i'),'.mat'],'acute','contactlengths','contactcounts','sqdis1','sqdis2');

end




pars=[500 1000 2000 4000 8000 10000 12000 14000 16000 32000]/1000;
N=length(pars);

for p = parnrs
    
parnr=parnrs(p);
load([savefolder,'acute-contact-sqdis',num2str(parnr, '%02i'),'.mat']);

% vals=contactcounts;
% figure(1)
% mpar=mean(vals,2);
% stdpar=std(vals');
% errorbar(1:N,mpar,stdpar,'*black');
% %boxplot(vals')
% xlabel('stiffness (kPa)');
% ylabel('contactcounts')
% set(gca,'XTick',1:N)
% xlim([0 N+1])
% ylim([0 Inf])
% set(gca, 'XTickLabel', pars)
% saveas(gca,[savefolder,'contactcounts',num2str(parnr, '%02i'),'.eps'],'eps')
% saveas(gca,[savefolder,'contactcounts',num2str(parnr, '%02i'),'.png'],'png')
% 
% 
% vals=sum(contactlengths,3)./contactcounts;
% vals(isnan(vals))=0;
% figure(2)
% mpar=mean(vals,2);
% stdpar=std(vals');
% errorbar(1:N,mpar,stdpar,'*black');
% %boxplot(vals')
% xlabel('stiffness (kPa)');
% ylabel('contact duration (MCS)')
% set(gca,'XTick',1:N)
% xlim([0 N+1])
% ylim([0 Inf])
% set(gca, 'XTickLabel', pars)
% saveas(gca,[savefolder,'contactdurations',num2str(parnr, '%02i'),'.eps'],'eps')
% saveas(gca,[savefolder,'contactdurations',num2str(parnr, '%02i'),'.png'],'png')
% 
% 
vals=acute;
figure(3)
mpar=mean(vals,2);
stdpar=std(vals');
errorbar(1:N,mpar,stdpar,'*black');
%boxplot(vals')
xlabel('stiffness (kPa)');
ylabel('fraction of obtuse angles')
set(gca,'XTick',1:N)
xlim([0 N+1])
ylim([0 Inf])
set(gca, 'XTickLabel', pars)
% saveas(gca,[savefolder,'acute',num2str(parnr, '%02i'),'.eps'],'eps')
% saveas(gca,[savefolder,'acute',num2str(parnr, '%02i'),'.png'],'png')



sqdisanalysis = zeros(size(sqdis1,1),2*size(sqdis1,2),size(sqdis1,3));
sqdisanalysis(:,1:size(sqdis1,2),:)=sqdis1;
sqdisanalysis(:,size(sqdis1,2)+1:end,:)=sqdis2;




vals=squeeze(sqdisanalysis(stiffness,:,:));
s=size(vals);
subt = 1:100:s(2);
vals=vals(:,subt);
figure(4)
mpar=mean(vals,1);
stdpar=std(vals);
errorbar(loopvecsqdis(subt),mpar,stdpar,'*blue');
%boxplot(vals')
xlim([loopvecsqdis(1) loopvecsqdis(end)])
ylim([0 Inf])
xlabel('Time (MCS)');
ylabel('MSD')
saveas(gca,[savefolder,'sqdis',num2str(parnr, '%02i'),'.eps'],'eps')
saveas(gca,[savefolder,'sqdis',num2str(parnr, '%02i'),'.png'],'png')


meansqdis=squeeze(mean(sqdisanalysis,2));
vals=pars;
N=length(vals);
x=loopvecsqdis;
M=meansqdis;
xl='time';
yl='ensemble MSD';
l='Youngs';
subt=1:100:5000;

figure(5)
ColorSet=varycolor(6);
set(gca,'ColorOrder',ColorSet)
hold all;
for i=1:N
    if(i<N/2+1)
plot(x(subt),M(i,subt),'LineWidth',2)
    else
plot(x(subt),M(i,subt),'--','LineWidth',2)
    end
end
xlabel(xl)
ylabel(yl)
h=legend('500' ,'1000' ,'2000', '4000', '8000', '10000', '12000', '14000', '16000', '32000','Location','NorthWest')
v = get(h,'title');
set(v,'string','stiffness (Pa)');
saveas(gca,[savefolder,'msd',num2str(parnr, '%02i'),'.eps'],'eps')
saveas(gca,[savefolder,'msd',num2str(parnr, '%02i'),'.png'],'png')
clf(gcf)

tstart=500;
tstart=find(x==tstart);
tend=5000;
tend=find(x==tend);
dispcoeffs=[];
confints=[];
for y = 1:N
    ft = fittype('poly1');    
    [fo, gof] = fit((0:tstart)', meansqdis(y,1:tstart+1)',ft, 'Lower',[0 0],'Upper',[100 100]);
    ci = confint(fo);
    dispcoeffs(y,1)=fo.p1/4;
    confints(y,1:2)=ci(:,1)'/4;
    
    [fo, gof,boe] = fit((tstart+1:tend)', meansqdis(y,tstart+1:tend)',ft,'Lower',[0 0],'Upper',[100 100])
    ci = confint(fo);
    dispcoeffs(y,2)=fo.p1/4;
    confints(y,3:4)=ci(:,1)'/4;
end


figure(6)
mpar=dispcoeffs(:,1);
L=confints(:,1)-dispcoeffs(:,1);
U=confints(:,2)-dispcoeffs(:,1);
errorbar(1:N,mpar,L,U,'*black');
xlabel('stiffness (kPa)');
ylabel('D')
set(gca,'XTick',1:N)
xlim([0 N+1])
set(gca, 'XTickLabel', pars)
saveas(gca,[savefolder,'Dfirst',num2str(parnr, '%02i'),'.eps'],'eps')
saveas(gca,[savefolder,'Dfirst',num2str(parnr, '%02i'),'.png'],'png')

figure(7)
mpar=dispcoeffs(:,2);
L=confints(:,3)-dispcoeffs(:,2);
U=confints(:,4)-dispcoeffs(:,2);
errorbar(1:N,mpar,L,U,'*black');
xlabel('stiffness (kPa)');
ylabel('D')
set(gca,'XTick',1:N)
xlim([0 N+1])
set(gca, 'XTickLabel', pars)
saveas(gca,[savefolder,'Dsecond',num2str(parnr, '%02i'),'.eps'],'eps')
saveas(gca,[savefolder,'Dsecond',num2str(parnr, '%02i'),'.png'],'png')



end



aantal=sum(acute(5,:)*1);
aantal=floor(aantal)
1-binocdf(aantal,100*1,0.25)

m=30;
r=randi(100,m,1);
format shortEng
format
ps1=[];
for i = 1:N
[h,ps1(i)]=ttest2(acute(4,:),acute(i,:),0.05,'both','equal');
end
ps1

ps1=[];
for i = 1:N
[h,ps1(i)]=ttest2(acute(4,r),acute(i,r),0.05,'both','unequal');
end
ps1

ps1=[];
for i = 1:N
[h,ps1(i)]=ztest(acute(4,:)-acute(i,:),0,std(acute(4,:),acute(i,:)))
end
ps1




ps2=[];
for i = 1:N
[h,ps2(i)]=ttest2(contactcounts(4,:),contactcounts(i,:));
end
ps2

[h,p]=kstest((acute(4,:)-mean(acute(4,:)))./std(acute(4,:)))
[h,p]=vartest2(acute(4,:),acute(5,:))
