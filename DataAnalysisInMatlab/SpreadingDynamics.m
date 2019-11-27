addpath('/path/to/functions/')

folderstart1 = '/path/to/Data_Pars';
Model='/path/to/model/';
checkexistfile='cpmfem00000.png';
savefolder='/path/to/savefolder/';
mkdir(savefolder)

par1=[1000,5000,10000,20000,50000,100000,100000000]; %youngs

par2 = [3 4 5 6 ]; %lambdafa
%par2 = [2500 2750 3000 3250 3500]; %jcm
%par2 = [4000000 6000000 8000000 10000000 20000000]; %capacityfa
par2=[1 1.5 2 2.5 3]; %temp
%par2 = [0.01 0.025 0.05 0.08 0.1]; %growthfa
par2 = [ 0 1 2 3]; %plaque
%par2 = [0.1e-9 0.25e-9 1e-9 5e-9 10e-9]*100; %v
%par2 =[0.00025 0.0005 0.001 0.0015 0.002]; %lrtension
%par2 = [2000 5000 10000 20000;] %confstress]

%par2=[1 1.5 2 2.5 3]; %temp


%par2=[100 500 1000 1500 2000] %pderepeat for tmcs




%par2=5; %lffa

%par2 =[0 1]; %im

sim=25;
M=size(par1,2);
N=size(par2,2);
MCS1=0;
MCS2=2000;
tstep=10; %dit moet 25 zijn
loopvec=MCS1:tstep:MCS2;
areas=zeros(length(loopvec),M,N,sim);
eccs=zeros(length(loopvec),M,N,sim);
rpas=zeros(length(loopvec),M,N,sim);




folder1=folderstart1
notexist=[];


for m = 1:M
    m
    for n = 1:N
        n
        for s = 1:sim
            s
            subfolder1 = [folder1,num2str(m,'%03i'),'-',num2str(n, '%03i'),'Sim',num2str(s, '%03i'),'/'];
            %if(~exist([subfolder1,checkexistfile]))
            %    cd(Model)
            %    system(['./CPMFEM ',subfolder1])
            %end
            if(exist([subfolder1,checkexistfile]))
            areas(:,m,n,s)=loadvariableMCS([subfolder1,'area.txt'],loopvec,1);
            eccs(:,m,n,s)=loadvariableMCS([subfolder1,'ecc.txt'],loopvec,1);
            %rpas(:,m,n,s)=loadvariableMCS([subfolder1,'ratiopa.txt'],loopvec,1);

            end
            if(~exist([subfolder1,checkexistfile]))
               notexist=[notexist;[m,n,s]]
            end
        end
    end
end
notexist

save([savefolder,'areas.mat'],'areas');
save([savefolder,'ecc.mat'],'eccs');


%load([savefolder,'areas.mat'])
%load([savefolder,'ecc.mat'])


areas=areas*6.25;
areas(areas==0)=NaN;

meanarea=nanmean(areas,4);
stdarea=nanstd(areas,0,4);
meanecc=nanmean(eccs,4);
stdecc=nanstd(eccs,0,4);
meanrpa=nanmean(rpas,4);
stdrpa=nanstd(rpas,0,4);
%errorbar(1:M,meanarea(end,:,3),stdarea(end,:,3),'*r','LineWidth',3);
%set(gca,'XTick',par1)
%set(gca,'FontSize',14)
%set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
%saveas(gca,[savefolder,'areas.eps'],'epsc')
%saveas(gca,[savefolder,'areas.png'],'png')



figure(1)
barwitherr(squeeze(stdarea(end,:,:)),squeeze(meanarea(end,:,:)))
set(gca,'XTickLabel',{'1','5','10','20','50','100','100000000'})
 %legendCell = cellstr(num2str(par2', '%-d'))
%hleg=legend(legendCell,'Location','NorthWest')
ylim([0 10000])
xlabel('Youngs modulus (kPa)')
% hlt = text(...
%     'Parent', hleg.DecorationContainer, ...
%     'String', 'values', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Position', [0.5, 1.05, 0], ...
%     'Units', 'normalized','FontSize',12,'fontWeight','bold')
ylabel('cell area (μm^2)')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure1.png'],'png')
saveas(gca,[savefolder,'figure1.eps'],'epsc')


figure(1)
pars=par1;
pars(7)=1.1*par1(6);
shadedErrorBar(pars/1000,squeeze(meanarea(end,:,1)),squeeze(stdarea(end,:,1)),'--o')
xticks(pars/1000)
set(gca,'XTickLabel',{'1','5','10','20','50','100','100000000'})
 %legendCell = cellstr(num2str(par2', '%-d'))
%hleg=legend(legendCell,'Location','NorthWest')
ylim([0 7500])
xlabel('Youngs modulus (kPa)')
% hlt = text(...
%     'Parent', hleg.DecorationContainer, ...
%     'String', 'values', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Position', [0.5, 1.05, 0], ...
%     'Units', 'normalized','FontSize',12,'fontWeight','bold')
ylabel('cell area (μm^2)')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure1-2.png'],'png')
saveas(gca,[savefolder,'figure1-2.eps'],'epsc')

figure(1)
colors=varycolor(N)
pars=par1;
pars(7)=1.1*par1(6);
for n=1:N
shadedErrorBar(pars/1000,squeeze(meanarea(end,:,n)),squeeze(stdarea(end,:,n)),{'--o','Color',colors(n,:),'LineWidth',2,'MarkerSize',10})
xticks(pars/1000)
hold on
end
set(gca,'XTickLabel',{'1','5','10','20','50','100','100000000'})
 %legendCell = cellstr(num2str(par2', '%-d'))
%hleg=legend(legendCell,'Location','NorthWest')
ylim([0 7500])
xlabel('Youngs modulus (kPa)')
% hlt = text(...
%     'Parent', hleg.DecorationContainer, ...
%     'String', 'values', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Position', [0.5, 1.05, 0], ...
%     'Units', 'normalized','FontSize',12,'fontWeight','bold')
ylabel('cell area (μm^2)')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure1-3.png'],'png')
saveas(gca,[savefolder,'figure1-3.eps'],'epsc')

figure(1)
colors=varycolor(N)
pars=par1;
pars(7)=1.1*par1(6);
for n=1:N
shadedErrorBar(pars/1000,squeeze(meanecc(end,:,n)),squeeze(stdecc(end,:,n)),{'--o','Color',colors(n,:),'LineWidth',2,'MarkerSize',10})
xticks(pars/1000)
hold on
end
set(gca,'XTickLabel',{'1','5','10','20','50','100','100000000'})
 %legendCell = cellstr(num2str(par2', '%-d'))
%hleg=legend(legendCell,'Location','NorthWest')
ylim([0 1.2])
xlabel('Youngs modulus (kPa)')
% hlt = text(...
%     'Parent', hleg.DecorationContainer, ...
%     'String', 'values', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Position', [0.5, 1.05, 0], ...
%     'Units', 'normalized','FontSize',12,'fontWeight','bold')
ylabel('cell eccentricity')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure2-3.png'],'png')
saveas(gca,[savefolder,'figure2-3.eps'],'epsc')





% figure(1)
% testarea=areas(51:201,:,:,:);
% testarea=squeeze(nanmean(testarea,1));
% stdareaMCS=nanstd(testarea,0,3);
% meanareaMCS=nanmean(testarea,3);
% barwitherr(stdareaMCS,meanareaMCS)
% set(gca,'XTickLabel',{'0.5','1','5','10','20','1000000'})
%  legendCell = cellstr(num2str(par2', '%-d'))
% hleg=legend(legendCell,'Location','NorthWest')
% ylim([0 7500])
% xlabel('Youngs modulus (kPa)')
% hlt = text(...
%     'Parent', hleg.DecorationContainer, ...
%     'String', 'values', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Position', [0.5, 1.05, 0], ...
%     'Units', 'normalized','FontSize',12,'fontWeight','bold')
% ylabel('cell area')
% set(gca,'FontSize',12)
% set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
% saveas(gca,[savefolder,'figure1.png'],'png')

figure(2)
barwitherr(squeeze(stdecc(end,:,:)),squeeze(meanecc(end,:,:)))
set(gca,'XTickLabel',{'1','5','10','20','50','100','100000000'})
 legendCell = cellstr(num2str(par2', '%-d'))
hleg=legend(legendCell,'Location','NorthWest')
ylim([0 1])
xlabel('Youngs modulus (kPa)')
hlt = text(...
    'Parent', hleg.DecorationContainer, ...
    'String', 'values', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Position', [0.5, 1.05, 0], ...
    'Units', 'normalized','FontSize',12,'fontWeight','bold')
ylabel('cell eccentricity')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure2.png'],'png')
saveas(gca,[savefolder,'figure2.eps'],'epsc')

figure(2)
pars=par1;
pars(7)=1.1*par1(6);
shadedErrorBar(pars/1000,squeeze(meanecc(end,:,2)),squeeze(stdecc(end,:,2)),'--o')
xticks(pars/1000)
set(gca,'XTickLabel',{'1','5','10','20','50','100','100000000'})
 %legendCell = cellstr(num2str(par2', '%-d'))
%hleg=legend(legendCell,'Location','NorthWest')
ylim([0 1.2])
xlabel('Youngs modulus (kPa)')
% hlt = text(...
%     'Parent', hleg.DecorationContainer, ...
%     'String', 'values', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Position', [0.5, 1.05, 0], ...
%     'Units', 'normalized','FontSize',12,'fontWeight','bold')
ylabel('ecc')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure2-1.png'],'png')
saveas(gca,[savefolder,'figure2-1.eps'],'epsc')



tend=find(loopvec==MCS2);

colors=varycolor(M);
figure(3)
n=4;
for m=1:M
    shadedErrorBar(loopvec(2:1:tend)*10/3600,meanarea(2:1:tend,m,n),stdarea(2:1:tend,m,n),{'Color',colors(m,:),'LineWidth',2,'MarkerSize',10})
    hold on
end
ylim([0 7500])
xlim([0 2000*10/60/60])
%legendCell = cellstr(num2str(par1', '%-d'))
%hleg=legend(legendCell,'Location','NorthWest')
xlabel('time (h)')
ylabel(['cell area (μm^2)'])    
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
saveas(gca,[savefolder,'figure3.png'],'png')
saveas(gca,[savefolder,'figure3.eps'],'epsc')


hold off
ft=fittype('A*(erf((t-t0)/tau)+1)', 'independent', 't', 'dependent', 'y' );
test=squeeze(meanarea(:,5,4));
fitresult=fit(loopvec',test,ft,'StartPoint',[3000 500 100])
plot(loopvec,test);
hold on
y=fitresult.A*(erf((loopvec-fitresult.t0)/fitresult.tau)+1);
plot(loopvec,y);
fitresult.t0*10/60 
'minutes'

% colors=varycolor(M);
% figure(11)
% n=2;
% for m=1:M
%     shadedErrorBar(loopvec(2:1:tend),meanrpa(2:1:tend,m,n),stdrpa(2:1:tend,m,n),{'Color',colors(m,:),'LineWidth',2,'MarkerSize',10})
%     hold on
% end
% xlabel('time (MCS)')
% ylabel(['rpa'])    
% set(gca,'FontSize',14)
% set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
% saveas(gca,[savefolder,'figure10.png'],'png')
% 
colors=varycolor(N);
figure(9)
m=6;
for n=1:N
    shadedErrorBar(loopvec(2:1:tend),meanecc(2:1:tend,m,n),stdecc(2:1:tend,m,n),{'Color',colors(n,:),'LineWidth',2,'MarkerSize',10})
    hold on
    pause(1)
end
xlabel('time (MCS)')
ylabel(['ecc'])    
%hleg=legend('Location','NorthWest')
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
saveas(gca,[savefolder,'figure9.png'],'png')

colors=varycolor(M);
figure(4)
n=2;
for m=1:M
    shadedErrorBar(loopvec(2:1:tend),meanecc(2:1:tend,m,n),stdecc(2:1:tend,m,n),{'Color',colors(m,:),'LineWidth',2,'MarkerSize',10})
    hold on
    pause(1)
end
xlabel('time (MCS)')
ylabel(['cell eccentricity']) 
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
saveas(gca,[savefolder,'figure4.png'],'png')


fa=[];
fa2=[];
for m = 1:M
    m
    for n = 1:N
        n
        fa3=[];
        for s = 1:sim
            s
            subfolder1 = [folder1,num2str(m,'%03i'),'-',num2str(n, '%03i'),'Sim',num2str(s, '%03i'),'/'];
            fa2=load([subfolder1,'fa/fa02000.txt']);
            fa3=[fa3;fa2(:)];
        end
       fa(:,m,n)=fa3(:);

    end
end
save([savefolder,'fa.mat'],'fa');

load([savefolder,'fa.mat'])

%fa=fa*8e-15/50*1e12;

figure(8)
testfa=squeeze(fa(:,:,1));
testfa(testfa<5000)=NaN;
testfa(testfa==0)=NaN;
% boxplot(testfa)
% set(gca,'XTickLabel',{'0.5','1','5','10','20','50','100','1000000'})
% ylabel('N')
% nrfa=sum(~isnan(testfa))
% saveas(gca,[savefolder,'figure8.png'],'png')
% 
% barwitherr(nanstd(testfa),nanmean(testfa))
% nanstd(testfa)./nanmean(testfa)

hist(testfa(:,1:7),20)
ylim([0 2500])
xlim([5000 13000])
xlabel('FA area (μm^2)')
ylabel('count')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
saveas(gca,[savefolder,'figure11-2.png'],'png')
saveas(gca,[savefolder,'figure11-2.eps'],'epsc')


[counts, bins] = hist(testfa(:,1:7),20);
plot(bins, counts)


A=[];
A.one=testfa(:,1);
A.two=testfa(:,2);
A.three=testfa(:,3);
A.four=testfa(:,4);
A.five=testfa(:,5);
A.six=testfa(:,6);
A.seven=testfa(:,7);

nhist(A,'median','noerror','smooth','numbers','samebins')
saveas(gca,[savefolder,'figure11-4.png'],'png')
saveas(gca,[savefolder,'figure11-4.eps'],'epsc')




