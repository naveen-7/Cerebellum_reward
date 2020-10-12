

%% GET ALL STATS FOR POPULATION
% written by naveen at cumc on 8/12/17

% CER_analysis_03g_LINEAR

load(POP_file)

% Start_T = -450;
% End_T  =1050;
% StartLIM_T = -400; EndLIM_T = 800;
% 
% Start_M = -750;
% End_M  =750;
% StartLIM_M = -500; EndLIM_M = 700;

%%% NEW  %%% 11/27/19
Start_T = -1550;
End_T  =1550;
StartLIM_T = -1400; EndLIM_T = 700;

Start_M = -1550;
End_M  =1550;
StartLIM_M = -700; EndLIM_M = 1400;


WN_T = W_PREV_T;
CN_T = C_PREV_T;

WN_M = W_PREV_M;
CN_M = C_PREV_M;

time_T = Start_T:End_T;


NW_M = W_CURR_M;
NC_M = C_CURR_M;
time_M = Start_M:End_M;








W_PREV_T_NORM = W_PREV_T-GLOBAL.MIN;
W_PREV_T_NORM = W_PREV_T_NORM/GLOBAL.MAX;

C_PREV_T_NORM = C_PREV_T-GLOBAL.MIN;
C_PREV_T_NORM = C_PREV_T_NORM/GLOBAL.MAX;

W_PREV_M_NORM = W_PREV_M-GLOBAL.MIN;
W_PREV_M_NORM = W_PREV_T_NORM/GLOBAL.MAX;

C_PREV_M_NORM = C_PREV_M-GLOBAL.MIN;
C_PREV_M_NORM = C_PREV_M_NORM/GLOBAL.MAX;


W_CURR_M_NORM = W_CURR_M-GLOBAL.MIN;
W_CURR_M_NORM = W_CURR_M_NORM/GLOBAL.MAX;

C_CURR_M_NORM = C_CURR_M-GLOBAL.MIN;
C_CURR_M_NORM = C_CURR_M_NORM/GLOBAL.MAX;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:4
    CW_DELTA.statsWN(index) = NaN;
    CW_DELTA.statsCN(index) = NaN;
    CW_DELTA.statsNW(index) = NaN;
    CW_DELTA.statsNC(index) = NaN;
    
    CW_DELTA.statsWC(index) = NaN;
    CW_DELTA.statsWW(index) = NaN;
    CW_DELTA.statsCW(index) = NaN;
    CW_DELTA.statsCC(index) = NaN; 
    
    CW_DELTA.statsMT(index) = NaN; 
    
    CW_DELTA.statsP(index) = NaN; 
     
end

for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        % avg ----------------------------------------------
        CW_DELTA.statsWN(index) = nanmean(W_PREV_T(1,INDS));
        CW_DELTA.statsCN(index) = nanmean(C_PREV_T(1,INDS));
        CW_DELTA.statsNW(index) = NaN;
        CW_DELTA.statsNC(index) = NaN;
        
        %         % peak ---------------------------------------------
        %         CW_DELTA.statsWN(index) = nanmax(W_PREV_T(1,INDS));
        %         CW_DELTA.statsCN(index) = nanmax(C_PREV_T(1,INDS));
        %         CW_DELTA.statsNW(index) = NaN;
        %         CW_DELTA.statsNC(index) = NaN;
        clear INDS
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        CW_DELTA.statsWC(index) = nanmean(WC_T(INDS));
        CW_DELTA.statsWW(index) = nanmean(WW_T(INDS));
        CW_DELTA.statsCW(index) = nanmean(CW_T(INDS));
        CW_DELTA.statsCC(index) = nanmean(CC_T(INDS));
        
        
        
        
        %%%%% norm stuff %%%%
        
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        % avg ----------------------------------------------
        CW_DELTA.statsWN(index) = nanmean(W_PREV_T(1,INDS));
        CW_DELTA.statsCN(index) = nanmean(C_PREV_T(1,INDS));
        CW_DELTA.statsNW(index) = NaN;
        CW_DELTA.statsNC(index) = NaN;
        
        
        
        
        clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsP(index) = ttest_NN(W_PREV_T_NORM(1,INDS),C_PREV_T_NORM(1,INDS));
       
        
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    clear INDS
    INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
    % avg ----------------------------------------------
    CW_DELTA.statsWN(index) = nanmean(W_PREV_M(1,INDS));
    CW_DELTA.statsCN(index) = nanmean(C_PREV_M(1,INDS));
    CW_DELTA.statsNW(index) = NaN;
    CW_DELTA.statsNC(index) = NaN;
    
    %     % peak ---------------------------------------------
    %     CW_DELTA.statsWN(index) = nanmax(W_PREV_T(1,INDS));
    %     CW_DELTA.statsCN(index) = nanmax(C_PREV_T(1,INDS));
    %     CW_DELTA.statsNW(index) = NaN;
    %     CW_DELTA.statsNC(index) = NaN;
    clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsWC(index) = nanmean(WC_M(INDS));
    CW_DELTA.statsWW(index) = nanmean(WW_M(INDS));
    CW_DELTA.statsCW(index) = nanmean(CW_M(INDS));
    CW_DELTA.statsCC(index) = nanmean(CC_M(INDS));
    
    
    
    clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsP(index) = ttest_NN(W_PREV_M_NORM(1,INDS),C_PREV_M_NORM(1,INDS));
    
    
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    clear INDS
    INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
    % avg ------------------------------------------
    CW_DELTA.statsWN(index) = NaN;
    CW_DELTA.statsCN(index) = NaN;
    CW_DELTA.statsNW(index) = nanmean(W_CURR_M(1,INDS));
    CW_DELTA.statsNC(index) = nanmean(C_CURR_M(1,INDS));
    
    %     % peak ----------------------------------------
    %     CW_DELTA.statsWN(index) = NaN;
    %     CW_DELTA.statsCN(index) = NaN;
    %     CW_DELTA.statsNW(index) = nanmax(W_CURR_M(1,INDS));
    %     CW_DELTA.statsNC(index) = nanmax(C_CURR_M(1,INDS));
    clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsWC(index) = nanmean(WC_M(INDS));
    CW_DELTA.statsWW(index) = nanmean(WW_M(INDS));
    CW_DELTA.statsCW(index) = nanmean(CW_M(INDS));
    CW_DELTA.statsCC(index) = nanmean(CC_M(INDS));
    
    
   
    clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsP(index) = ttest_NN(W_CURR_M_NORM(1,INDS),C_CURR_M_NORM(1,INDS));
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:4
    CW_DELTA.statsMT(index) = NaN;
end

for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        CW_DELTA.statsMT(index) = nanmean(MT_T(1,INDS));
    end
end

for index=3:4
    if DELTA.EpochFlag(index)==1
        Time = Start_M:End_M;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        CW_DELTA.statsMT(index) = nanmean(MT_M(1,INDS));
    end
end

save(POP_file,'CW_DELTA','-append');
save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');




%%%%%%%%%%%%%% all stddev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:4
    CW_DELTA.statsWNstd(index) = NaN;
    CW_DELTA.statsCNstd(index) = NaN;
    CW_DELTA.statsNWstd(index) = NaN;
    CW_DELTA.statsNCstd(index) = NaN;
    
    CW_DELTA.statsWCstd(index) = NaN;
    CW_DELTA.statsWWstd(index) = NaN;
    CW_DELTA.statsCWstd(index) = NaN;
    CW_DELTA.statsCCstd(index) = NaN; 
    
    CW_DELTA.statsMTstd(index) = NaN; 
    
    CW_DELTA.statsPstd(index) = NaN; 
     
end

for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        % avg ----------------------------------------------
        CW_DELTA.statsWNstd(index) = nanstd(W_PREV_T(1,INDS));
        CW_DELTA.statsCNstd(index) = nanstd(C_PREV_T(1,INDS));
        CW_DELTA.statsNWstd(index) = NaN;
        CW_DELTA.statsNCstd(index) = NaN;
        
        %         % peak ---------------------------------------------
        %         CW_DELTA.statsWN(index) = nanmax(W_PREV_T(1,INDS));
        %         CW_DELTA.statsCN(index) = nanmax(C_PREV_T(1,INDS));
        %         CW_DELTA.statsNW(index) = NaN;
        %         CW_DELTA.statsNC(index) = NaN;
        clear INDS
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        CW_DELTA.statsWCstd(index) = nanstd(WC_T(INDS));
        CW_DELTA.statsWWstd(index) = nanstd(WW_T(INDS));
        CW_DELTA.statsCWstd(index) = nanstd(CW_T(INDS));
        CW_DELTA.statsCCstd(index) = nanstd(CC_T(INDS));
        
        
        
        
        %%%%% norm stuff %%%%
        
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        % avg ----------------------------------------------
        CW_DELTA.statsWNstd(index) = nanstd(W_PREV_T(1,INDS));
        CW_DELTA.statsCNstd(index) = nanstd(C_PREV_T(1,INDS));
        CW_DELTA.statsNWstd(index) = NaN;
        CW_DELTA.statsNCstd(index) = NaN;
        
        
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    clear INDS
    INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
    % avg ----------------------------------------------
    CW_DELTA.statsWNstd(index) = nanstd(W_PREV_M(1,INDS));
    CW_DELTA.statsCNstd(index) = nanstd(C_PREV_M(1,INDS));
    CW_DELTA.statsNWstd(index) = NaN;
    CW_DELTA.statsNCstd(index) = NaN;
    
    %     % peak ---------------------------------------------
    %     CW_DELTA.statsWN(index) = nanmax(W_PREV_T(1,INDS));
    %     CW_DELTA.statsCN(index) = nanmax(C_PREV_T(1,INDS));
    %     CW_DELTA.statsNW(index) = NaN;
    %     CW_DELTA.statsNC(index) = NaN;
    clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsWCstd(index) = nanstd(WC_M(INDS));
    CW_DELTA.statsWWstd(index) = nanstd(WW_M(INDS));
    CW_DELTA.statsCWstd(index) = nanstd(CW_M(INDS));
    CW_DELTA.statsCCstd(index) = nanstd(CC_M(INDS));
   
    
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    clear INDS
    INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
    % avg ------------------------------------------
    CW_DELTA.statsWNstd(index) = NaN;
    CW_DELTA.statsCNstd(index) = NaN;
    CW_DELTA.statsNWstd(index) = nanstd(W_CURR_M(1,INDS));
    CW_DELTA.statsNCstd(index) = nanstd(C_CURR_M(1,INDS));
    
    %     % peak ----------------------------------------
    %     CW_DELTA.statsWN(index) = NaN;
    %     CW_DELTA.statsCN(index) = NaN;
    %     CW_DELTA.statsNW(index) = nanmax(W_CURR_M(1,INDS));
    %     CW_DELTA.statsNC(index) = nanmax(C_CURR_M(1,INDS));
    clear INDS
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.statsWCstd(index) = nanstd(WC_M(INDS));
    CW_DELTA.statsWWstd(index) = nanstd(WW_M(INDS));
    CW_DELTA.statsCWstd(index) = nanstd(CW_M(INDS));
    CW_DELTA.statsCCstd(index) = nanstd(CC_M(INDS));
    
 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:4
    CW_DELTA.statsMTstd(index) = NaN;
end

for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        CW_DELTA.statsMTstd(index) = nanstd(MT_T(1,INDS));
    end
end

for index=3:4
    if DELTA.EpochFlag(index)==1
        Time = Start_M:End_M;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        CW_DELTA.statsMTstd(index) = nanstd(MT_M(1,INDS));
    end
end


save(POP_file,'CW_DELTA','-append');
save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DELTA.EpochFlag(1)==0 & DELTA.EpochFlag(2)==0 & DELTA.EpochFlag(3)==0
    CW_DELTA.PREV=[];
    DELTA.HELPER(1)=NaN;
end

if DELTA.EpochFlag(4)==0
    CW_DELTA.CURR=[];
    DELTA.HELPER(2)=NaN;
end








save(POP_file,'CW_DELTA','-append');
save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');







% end
