savefolder = 'path/to/savefolder/TwoCells_J';
parnr=4

load([savefolder,'acute-contact-sqdis',num2str(parnr, '%02i'),'.mat']);
sqdisanalysis = zeros(size(sqdis1,1),2*size(sqdis1,2),size(sqdis1,3));
sqdisanalysis(:,1:size(sqdis1,2),:)=sqdis1;
sqdisanalysis(:,size(sqdis1,2)+1:end,:)=sqdis2;

savefolder = 'path/to/savefolder/SingleCell_Lambda/';
load([savefolder,'length-area-ecc-sqdis',num2str(parnr, '%02i'),'.mat']);

savefolder = '/ufs/rens/CPM_FEM_Project/Results/ForPaperRene/TwoCells_J/';

sqdisend1c = squeeze(sqdis(:,:,end));
sqdisend2c = squeeze(sqdisanalysis(:,:,end));

ttestresult=[];
pvals =[];

for st = 1:10
    [ttestresult(st),pvals(st)]=ttest2(sqdisend1c(st,:),sqdisend2c(st,:));
end

for st = 1:10
    [ttestresult(st),pvals(st)]=ttest2(sqdisend1c(st,:),sqdisend2c(st,101:200));
end