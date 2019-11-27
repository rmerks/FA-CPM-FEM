function s = fstdangles(a,m,dim)
        ar=ones(1,length(size(a))-1);
        ar=[ar,size(a,dim)];
        M=repmat(m,ar);
        order1=1:length(size(a));
        order=order1;
        order(dim)=length(size(a));
        for k = 1:length(size(a))-dim;
            order(dim+k)=order1(dim+k-1);
        end
        M=permute(M,order);
        
        z1 = cos(a*pi/180);
        z2 = sin(a*pi/180);
        w1 = cos(M*pi/180);
        w2 = sin(M*pi/180);
        ab = acos(w1.*z1+w2.*z2)*180/pi;
    ab(ab>=90)=180-ab(ab>=90);
    S2 = sum(ab.^2,dim)/(size(a,dim)-1);
    s=squeeze(sqrt(S2));
    %using http://magician.ucsd.edu/essentials/WebBookse69.html 
    %estimated standard deviation
end