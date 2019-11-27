addpath('path/to/functions')
folderstart = '/path/to/folder/Data_Pars';
Model='/path/to/model';
checkexistfile='cpmfem00000.png'

savefolder = '/path/to/savefolder/';
mkdir(savefolder)

parnrs=1:4;


pars=10000:1000:16000;
pars=[0.5 1 2 4 10 12 14 16 32]*1000;
N=length(pars);
sim=100;


MCS1=0;
MCS2=5000;
tstep=1;
loopvec=MCS1:tstep:MCS2;

for p = 1:length(parnrs)
parnr=parnrs(p);

folder = [folderstart,num2str(parnr, '%03i'),'-'];


clength=zeros(N,sim,length(loopvec));


for n = 1:N
    for s = 1:sim
        s
        folder1 = [folder,num2str(n, '%03i'),'Sim',num2str(s, '%03i'),'/'];
        if(~exist([folder1,checkexistfile]))
            cd(Model)
            system(['./CPMFEM ',folder1])
        end
        clength(n,s,:)=loadvariableMCS([folder1,'length.txt'],loopvec,1)';
    end
end

save([savefolder,'clengths',num2str(parnr, '%02i'),'.mat'],'clength','pars');
end


for p = 1:length(parnrs)
parnr=parnrs(p);
load([savefolder,'clengthsT',num2str(parnr, '%02i'),'.mat']);
vals=squeeze(mean(clength(:,:,21:101),3));
figure(1)
mpar=mean(vals,2);
stdpar=std(vals');
errorbar(1:N,mpar,stdpar,'*black');
%boxplot(vals')
xlabel('stiffness (kPa)');
ylabel('cell length (pixels)')
set(gca,'XTick',1:N)
xlim([0 N+1])
set(gca, 'XTickLabel', pars/1000)
saveas(gca,[savefolder,'clength',num2str(parnr, '%02i'),'.eps'],'eps')
saveas(gca,[savefolder,'clength',num2str(parnr, '%02i'),'.png'],'png')
end
