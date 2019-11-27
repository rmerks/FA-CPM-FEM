function acute = GetObtuseAcute(folder,loopvec)

variable='angle';
file = [folder,variable,'.txt'];

M=loadvariableMCS(file,loopvec,2);

yesacute=[];
for m = 1:length(loopvec);
    [a(1),a(2),a(3)]=AngleTwoCell(M(m,:));
    if(a(3)>90)
    yesacute(m)=1;   
    end
end
acute = sum(yesacute)/length(yesacute);

end