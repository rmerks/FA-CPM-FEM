function allval = get_timeval(folder,variable,t)
folder=[folder,variable,'/'];
if(exist(folder))
    allval = load([folder,variable,num2str(t, '%05i'),'.out']);
else
    'Does not exist'
    allval=2^32;
end

end