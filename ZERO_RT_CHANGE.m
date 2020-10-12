

%% function ZERO_RT_CHANGE



RT_now = Infos(:,14);

RT_now_BEF=(RT_now(CHANGE-15:CHANGE-1));
RT_now_AFT=(RT_now(CHANGE:CHANGE+19));

if ttest_NN(RT_now_BEF,RT_now_AFT)<0.05
    ZERO_RT_FLAG=0;   % SIGNIFICANT 
end

if ttest_NN(RT_now_BEF,RT_now_AFT)>0.05
    ZERO_RT_FLAG=1   % NOT SIGNIFICANT 
end



corr = Infos(:,10);
corr(corr==2)=0;

LC_CHANGE(1) = nanmean(corr(CHANGE-15:CHANGE-1));
LC_CHANGE(2) = nanmean(corr(CHANGE:CHANGE+19));

RT_CHANGE(1) = nanmean(RT_now(CHANGE-15:CHANGE-1));
RT_CHANGE(2) = nanmean(RT_now(CHANGE:CHANGE+19));







save(POP_file,'ZERO_RT_FLAG','LC_CHANGE','RT_CHANGE','-append');
save(MERGE_file,'ZERO_RT_FLAG','LC_CHANGE','RT_CHANGE','-append');
save(ALLCELLS_file,'ZERO_RT_FLAG','LC_CHANGE','RT_CHANGE','-append');




%% end