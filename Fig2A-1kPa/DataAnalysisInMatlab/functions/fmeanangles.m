function m = fmeanangles(a,dim,half)
if(half)
cosa = cos(2*a*pi/180); %deal with half circular data, by multiplying by two
sina = sin(2*a*pi/180);
end
if(~half)
cosa = cos(a*pi/180); 
sina = sin(a*pi/180);
end
v1 = nanmean(cosa,dim);
v2 = nanmean(sina,dim);
m = atan2(v2,v1)*180/pi;
if(half)
m=m/2;
end
m = mod(m,180);
m=squeeze(m);
end