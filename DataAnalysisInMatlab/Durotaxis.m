addpath('/path/to/functions/')
folderstart1 = '/path/to/dur/Data_Pars';
Model='/path/to/model/';
checkexistfile='cpmfem00000.png';
savefolder='/path/to/savefolder/';
mkdir(savefolder)


%par1=0.0002;
par2=[10 20 50 80 100 200]; %gradient
%par2=0.0002;


par2=[10 20 50 80 100 200]; %gradient
par1=[6000];
par2=[1 2 3]; %temp
par1=[0.00015 0.0002 0.00025]; %lambda-area

%par1=1
%par2=1

sim=25;
M=size(par1,2);
N=size(par2,2);
MCS1=0;
MCS2=10000;
tstep=10; 
loopvec=MCS1:tstep:MCS2;

folder1=folderstart1;

csx=zeros(M,N,sim,length(loopvec));
csy=zeros(M,N,sim,length(loopvec));

notexist=[];

for m = 1:M
    m
    for n = 1:N
        n
        for s = 1:sim
            s
            tic
           subfolder1 = [folder1,num2str(m,'%03i'),'-',num2str(n, '%03i'),'Sim',num2str(s, '%03i'),'/'];
            if(~exist([subfolder1,checkexistfile]))
               notexist=[notexist;[m,n,s]]
            end
            for tt = 1:1:length(loopvec)
            mcs=loopvec(tt);
            if(exist([subfolder1,checkexistfile]))
             sigma=load([subfolder1,'sigmas/sigma',num2str(mcs, '%05i'),'.txt']);
             A=sigma;
             C=cellfun(@(n) 1:n, num2cell(size(A)),'uniformoutput',0);
                 [C{:}]=ndgrid(C{:});
                 C=cellfun(@(x) x(:), C,'uniformoutput',0);
                 C=[C{:}];
                 com2=A(:).'*C/sum(A(:),'double');
             csx(m,n,s,tt)=com2(1);
             csy(m,n,s,tt)=com2(2);
            end



            end
            toc
        end
    end
end
notexist

save([savefolder,'csx.mat'],'csx');
save([savefolder,'csy.mat'],'csy');


% load([savefolder,'csx.mat'])
% load([savefolder,'csy.mat'])


%um
csx=csx*2.5;
csy=csy*2.5;

csx(csx==0)=NaN;
%m=1,n=3 is default
mcsx=squeeze(nanmean(csx,3));
stdcsx=squeeze(nanstd(csx,0,3));
%mcsx0=squeeze(nanmean(csx0,3));
%stdcsx0=squeeze(nanstd(csx0,0,3));


m=1;
n=3;
colors=varycolor(2);
tstep=20;
shadedErrorBar(loopvec(1:tstep:length(loopvec))*10/3600,squeeze(mcsx(m,n,1:tstep:length(loopvec))),squeeze(stdcsx(m,n,1:tstep:length(loopvec))),{'Color',colors(1,:),'LineWidth',2,'MarkerSize',10});
%hold on
%shadedErrorBar(loopvec(1:tstep:length(loopvec))*10/3600,mcsx0(1:tstep:length(loopvec)),stdcsx0(1:tstep:length(loopvec)),{'Color',colors(2,:),'LineWidth',2,'MarkerSize',10})
saveas(gca,[savefolder,'figure1.eps'],'epsc')


barwitherr(squeeze(stdcsx(:,:,end)'),squeeze(mcsx(:,:,end))');
saveas(gca,[savefolder,'figure2.eps'],'epsc')

barwitherr(squeeze(stdcsx(:,:,end)),squeeze(mcsx(:,:,end)));
saveas(gca,[savefolder,'figure3.eps'],'epsc')

m=[squeeze(mcsx(1,3,end)),mcsx0(end)];
s=[squeeze(stdcsx(1,3,end)),stdcsx0(end)];
barwitherr(s,m)
saveas(gca,[savefolder,'figure4.eps'],'epsc')



m=1;
n=3;
%trajectories
csx1=squeeze(csx(m,n,:,:));
csy1=squeeze(csy(m,n,:,:));

p=1:50:length(loopvec);
S=10;
colors=varycolor(S);
for m=1:S
plot(csx1(m,p)',csy1(m,p)','LineWidth',2,'Color',colors(m,:))
hold on
plot(csx1(m,end)',csy1(m,end)','.black','MarkerSize',25)
end
%hold on
%plot(100,100,'.black','MarkerSize',25)
%hold on
%plot(100,100,'+black','MarkerSize',25)
%xlim([80 180]*2.5)
%ylim([80 180]*2.5)
xlim([100 200]*2.5)
ylim([0 200]*2.5)
saveas(gca,[savefolder,'figure3.eps'],'epsc')








m=1;
n=3;
colors=varycolor(6);
tstep=20;
for n =3:3
shadedErrorBar(loopvec(1:tstep:length(loopvec))*10/3600,squeeze(mcsx(n,1:tstep:length(loopvec))),squeeze(stdcsx(n,1:tstep:length(loopvec))),{'Color',colors(n,:),'LineWidth',2,'MarkerSize',10});
hold on
end
saveas(gca,[savefolder,'figure4.eps'],'epsc')





slope=diff(mcsx')';
mean(slope')
std(slope')

time=loopvec*10/3600;
slopes=[];
for n=1:6
c=polyfit(time,mcsx(n,:),1);
slopes(n)=c(1);
end

figure(5)
plot(par2/2.5,slopes,'-*')
saveas(gca,[savefolder,'figure5.eps'],'epsc')



%dur
slopes=[]
for m=1:3
    for n=1:3
    c=polyfit(time,squeeze(mcsx(m,n,:))',1);
    slopes(m,n)=c(1);
    end
end



%dur0
c=polyfit(time,mcsx',1);




time=loopvec*10/3600;
slopes=[];
for n=1:6
    for s=1:sim
c=polyfit(time,squeeze(csx(1,n,s,:))',1);
slopes(n,s)=c(1);
    end
end
mslope=mean(slopes,2);
stdslope=std(slopes,1,2);
figure(6)
shadedErrorBar(par2/2.5,mslope,stdslope,'-*')
saveas(gca,[savefolder,'figure6.eps'],'epsc')



time=loopvec*10/3600;
slopes=[];
for n=1:3
    for s=1:sim
        for m =1:3
c=polyfit(time,squeeze(csx(m,n,s,:))',1);
slopes(m,n,s)=c(1);
        end
    end
end
mslope=mean(slopes,3);
stdslope=std(slopes,1,3);

barwitherr([33.0143,33.8960],[370.2335,366.2027]);
