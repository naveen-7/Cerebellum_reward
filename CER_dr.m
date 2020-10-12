
clear Activity_m
count=0;
for i=CHANGE-19:CHANGE+20
    count=count+1;
    Activity_m{count,1} = Spikes.S{i,1}-Infos(i,11);
end

Outcome = Infos(CHANGE-19:CHANGE+20,10);

save(POP_file,'Activity_m','Outcome','-append');
save(MERGE_file,'Activity_m','Outcome','-append');
save(ALLCELLS_file,'Activity_m','Outcome','-append');

