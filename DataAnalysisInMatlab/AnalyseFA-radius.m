folder = 'path/to/pl/Data_Pars';
savefolder='/path/to/savefolder/';

MCS1=0;
MCS2=2000;
tstep=10;
loopvec=MCS1:tstep:MCS2;

fa=sparse(25*NVX*NVY,7);
dis=sparse(25*NVX*NVY,7);

NVX=200;
NVY=200;

    X=repmat(1:NVX,NVX,1)';
    Y=repmat(1:NVY,NVY,1);
    
tt=length(loopvec);
folder1=folder;

                    for s=1:25
                        for m=1:7
                mcs=loopvec(tt);
                mcs
                s
               
                subfolder1 = [folder1,num2str(m,'%03i'),'-',num2str(1, '%03i'),'Sim',num2str(s, '%03i'),'/'];   
                if(exist(subfolder1))
                fanew=reshape(load([subfolder1,'fa/fa',num2str(mcs, '%05i'),'.txt']),NVX*NVY,1);
                fa(s*NVX*NVY+1:(s+1)*NVX*NVY,m)=fanew;


                s1=load([subfolder1,'sigmas/sigma',num2str(mcs, '%05i'),'.txt'])/20000;

                
                A=s1;            
                C=cellfun(@(n) 1:n, num2cell(size(A)),'uniformoutput',0);
                [C{:}]=ndgrid(C{:});
                C=cellfun(@(x) x(:), C,'uniformoutput',0);
                C=[C{:}];
                com1=A(:).'*C/sum(A(:),'double');
                
                
                disnew=reshape(sqrt((X-com1(1)).^2+(Y-com1(2)).^2),NVX*NVY,1);
                dis(s*NVX*NVY+1:(s+1)*NVX*NVY,m)=disnew;

                
                
                p1=find(fa<5000);
                fa(p1)=0;
                dis(p1)=0;
                end
             

                        end
                    end

for m =1:7
scatter(dis(:,m),fa(:,m))
hold on
end
xlim([0 30])


gm=[];
gs=[];
maxdis=25;
h=1/2.5;
for k=1:maxdis/h;
    for m=1:7
        p=find(dis(:,m)>h*k & dis(:,m)<h*(k+1));

g=fa(p,m);

gm(m,k)=nanmean(g);
gs(m,k)=nanstd(g);
    end
end

addpath('/path/to/functions/')

colors=varycolor(7);
for i=1:7
shadedErrorBar((h:h:maxdis)*2.5,gm(i,:),gs(i,:),{'Color',colors(i,:),'LineWidth',2,'MarkerSize',10})
hold on
end



saveas(gca,[savefolder,'figure-radiusfa.png'],'png')
saveas(gca,[savefolder,'figure-radiusfa.eps'],'epsc')


%saveas(gca,[savefolder,'figure-radiusfa-1.png'],'png')
%saveas(gca,[savefolder,'figure-radiusfa-1.eps'],'epsc')





