
clear Activity_m
count=0;
for i=CHANGE-14:CHANGE+15
    count=count+1;
    Activity_m{count,1} = Spikes.S{i,1}-Infos(i,11);
end

Outcome = Infos(CHANGE-14:CHANGE+15,10);

save(POP_file,'Activity_m','Outcome','-append');
save(MERGE_file,'Activity_m','Outcome','-append');
save(ALLCELLS_file,'Activity_m','Outcome','-append');

