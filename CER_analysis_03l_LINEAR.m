

DELTA.LEN = DELTA.END-DELTA.START;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PSTH_T_all PSTH_M_all
time_t = -1500:1500;
time_m = -1500:1500;
for i=1:size(Infos,1)
    PSTH_T_all(i,:) = PSTH_ONE_n(Spikes.S{i,1},Infos(i,4),time_t(1),time_t(end),40);
    PSTH_M_all(i,:) = PSTH_ONE_n(Spikes.S{i,1},Infos(i,11),time_m(1),time_m(end),40);
end


CW_DELTA.RAW = NaN(4,4);
CW_DELTA.USE = NaN(2,4);
    


i = find(DELTA.EpochFlag==1)
figure();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get trial outcomes
for j=1:length(i)
    
    
    if i(j)<=3
        corr_all = find(Infos(CHANGE:CHANGE+30,10)==1)+CHANGE-1 + 1;
        wrng_all = find(Infos(CHANGE:CHANGE+30,10)==2)+CHANGE-1 + 1;
    end
    
    if i(j)==4
        corr_all = find(Infos(CHANGE:CHANGE+30,10)==1)+CHANGE-1;
        wrng_all = find(Infos(CHANGE:CHANGE+30,10)==2)+CHANGE-1;
    end
    
    Corr_first5 = corr_all(1:5); Wrng_first5 = wrng_all(1:5);
    Corr_last5 = corr_all(end-5:end); Wrng_last5 = wrng_all(end-5:end);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make PSTHs
    
    if i(j)<=2
        PSTH_corr_first5 =  nanmean(PSTH_T_all(Corr_first5,:));
        PSTH_corr_last5 =  nanmean(PSTH_T_all(Corr_last5,:));
        PSTH_wrng_first5 =  nanmean(PSTH_T_all(Wrng_first5,:));
        PSTH_wrng_last5 =  nanmean(PSTH_T_all(Wrng_last5,:));
        time = time_t;
    end
    
    if i(j)>=3
        PSTH_corr_first5 =  nanmean(PSTH_M_all(Corr_first5,:));
        PSTH_corr_last5 =  nanmean(PSTH_M_all(Corr_last5,:));
        PSTH_wrng_first5 =  nanmean(PSTH_M_all(Wrng_first5,:));
        PSTH_wrng_last5 =  nanmean(PSTH_M_all(Wrng_last5,:));
        time = time_m;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get values in delta
    
    % RAW:
    %          F     T      M      P
    % corr1
    % wrng1
    % corr2
    % wrng2
    
    
    % USE:
    %               F     T      M      P
    % wrng1-corr1
    % wrng2-corr2
    
    
    
    
    
    TIME = find(time==DELTA.START(i(j))):find(time==DELTA.END(i(j)));
    CW_DELTA.RAW(1,i(j)) = nanmean(PSTH_corr_first5(TIME));
    CW_DELTA.RAW(2,i(j)) = nanmean(PSTH_wrng_first5(TIME));
    CW_DELTA.RAW(3,i(j)) = nanmean(PSTH_corr_last5(TIME));
    CW_DELTA.RAW(4,i(j)) = nanmean(PSTH_wrng_last5(TIME));
    
    CW_DELTA.USE(1,i(j)) = CW_DELTA.RAW(2,i(j))-CW_DELTA.RAW(1,i(j)); % wrong-corr
    CW_DELTA.USE(2,i(j)) = CW_DELTA.RAW(4,i(j))-CW_DELTA.RAW(3,i(j)); % wrong-corr
    
    
    
    
    
    
    subplot(2,2,i(j))
    hold on;
    plot(PSTH_corr_first5,'-b');
    plot(PSTH_corr_last5,'--b');
    plot(PSTH_wrng_first5,'-r');
    plot(PSTH_wrng_last5,'--r');
    
    
    
end



save(POP_file,'CW_DELTA','DELTA','-append');
save(MERGE_file,'CW_DELTA','DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','DELTA','-append');

 
    
    
