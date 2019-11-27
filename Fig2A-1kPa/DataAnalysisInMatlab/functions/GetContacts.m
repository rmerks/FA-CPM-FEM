function [contactcount contactlengths] = GetContacts(folder,loopvec)
variable='twocellcontact';
file = [folder,variable,'.txt'];

C=loadvariableMCS(file,loopvec,1);

if(all(C)==0)
    contactcount=0;
    contactlenghts=[];
end
CC=bwconncomp(C);
contactcount=CC.NumObjects;
contactlengths =[];
for i = 1:CC.NumObjects
    contactlengths = [contactlengths,length(CC.PixelIdxList{i})];
end

if(all(C==0))
    contactcount=0;
    clearvars contactlengths;
    contactlengths=0;
end


end

