function dircomp = DirectionComponents(sigma,N)
dircomp=[];
Ixx=0;
Iyy=0;
Ixy=0;

if(~isempty(N))
k=0;
cmasses = [];
for n = N
    k=k+1;
    [ix iy]=find(sigma==n);
    cmasses(k,:)=[mean(ix),mean(iy)];
end

k=0;
for n=N
    k=k+1;
    for x=1:size(sigma,1)
        for y=1:size(sigma,2)
            if(sigma(x,y)==n)
                Ixx=Ixx+(cmasses(k,1)-x)^2;
                Iyy=Iyy+(cmasses(k,2)-y)^2;
                Ixy=Ixy-(cmasses(k,1)-x)*(cmasses(k,2)-y);
            end
        end
    end
    I=[Ixx,Ixy;Ixy,Iyy];
    [V D]=eigs(I);
    d=90-atan2(V(2,1),V(1,1))*180/pi;
    dircomp(k,1)=mod(d,180);
    
end
end
end
