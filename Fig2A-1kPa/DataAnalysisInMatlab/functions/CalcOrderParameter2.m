function S = CalcOrderParameter(sigma,r,angles)
%2D orientational order parameter

%s=6;
%sigma = 1:s*s;
%sigma=reshape(sigma,s,s);
%r=10;
%angles=rand(s*s,1)*2*pi; %random.... should be around zero: but is actually more close to 0.25, just like in magriets paper
%angles=[ones(1,s*s/2)*pi,ones(1,s*s/2)*pi]'+pi/1000*randn(size(angles)); %aligned, should be one: correct

celldirs=[cos(2*angles*pi/180),sin(2*angles*pi/180)]; %to deal with halve circular data

x=1:size(sigma,1);
y=1:size(sigma,2);
ncells = length(angles);

%disp('calculate center of masses')

cmasses = [];
for n = 1:ncells
    [ix iy]=find(sigma==n);
    cmasses(n,:)=[mean(ix),mean(iy)];
end



%disp('calculate order parameter')

thetas=[];
for n = 1:ncells
    meandir=[0 0]; %what to do when no surrounding cells?
    lmeandir=1;
    X=repmat(x,size(sigma,2),1)';
    Y=repmat(y,size(sigma,1),1);
    disx=cmasses(n,1)-X;
    disy=cmasses(n,2)-Y;
    ids = find(disx.^2+disy.^2 <= r^2);
    boe = find(sigma(ids)>0); 
    if(~isempty(boe))
        ba = ids(boe);
        cells=sigma(ba);
        %length(unique(cells))
        %pause(0.5)
        meandir=mean(celldirs(cells',:),1);
        lmeandir=norm(meandir,2);
    end
%0.5*atan2(meandir(2),meandir(1))*180/pi
    if(~all(meandir==0))
    thetas=[thetas,acos(dot(meandir,celldirs(n,:))/lmeandir)*180/pi];
    end
    
end

thetas=thetas/2; %back to normal

%S=mean(3/2*cos(thetas*pi/180).^2-1/2); %3D-oops
S=mean(cos(2*thetas*pi/180)); %2D




end
