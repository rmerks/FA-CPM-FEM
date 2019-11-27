function M = loadvariableMCS(file,MCS,NumCols)
fileID=fopen(file);
FormatString=repmat('%f ',1,NumCols+1);
M=textscan(fileID,FormatString,'Headerlines',0);
M=cell2mat(M);
columns=[];
for m = MCS
    cid = find(M(:,1)==m);
    columns = [columns;cid];
end
M=M(columns,2:end);
fclose(fileID);
end