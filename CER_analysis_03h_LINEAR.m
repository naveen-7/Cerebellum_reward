
%% function CER_analysis_03h_LINEAR
% doing the memory analysis for each hand ------------------

% written by naveen at cumc on 8/12/17
%







SIGNAL = Spikes.S;


WRG_COLOUR = [255 128 128]/255;
COR_COLOUR = [0 179 179]/255;
OT_COLOUR = [81 81 81]/255;


ISWRNG_PREV = NaN(length(Infos),1);
for i=2:length(Infos)
    if (Infos(i-1,10)==2)
        ISWRNG_PREV(i,1)=1;
    end
end

ISCORR_PREV = NaN(length(Infos),1);
for i=2:length(Infos)
    if (Infos(i-1,10)==1)
        ISCORR_PREV(i,1)=1;
    end
end

ISWRNG_CURR = NaN(length(Infos),1);
for i=1:length(Infos)
    if (Infos(i,10)==2)
        ISWRNG_CURR(i,1)=1;
    end
end

ISCORR_CURR = NaN(length(Infos),1);
for i=1:length(Infos)
    if (Infos(i,10)==1)
        ISCORR_CURR(i,1)=1;
    end
end


Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;




IND_L = find(WHICHHAND==1);
IND_R = find(WHICHHAND==2);



for ind_hand=1:2
    
    if  ind_hand==1
        IND_H=IND_L;
        HAND = 'left';
    end
    
    if  ind_hand==2
        IND_H=IND_R;
        HAND = 'right';
    end
    
    
    FF = figure;
    SIG =40;
    
    
    % Prev Wrng
    IND = intersect(find(ISWRNG_PREV(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_H);
    W_PREV_NUM=length(IND);
    subplot(2,2,1)
    W_PREV_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
    subplot(2,2,2)
    W_PREV_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim;
    
    
    % Prev Corr
    IND = intersect(find(ISCORR_PREV(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_H);
    C_PREV_NUM=length(IND);
    subplot(2,2,1)
    C_PREV_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim; ylabel('Prev','fontsize',12);
    subplot(2,2,2)
    C_PREV_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;
    
    
    % Curr Wrng
    IND = intersect(find(ISWRNG_CURR(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_H);
    W_CURR_NUM=length(IND);
    subplot(2,2,3)
    W_CURR_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
    subplot(2,2,4)
    W_CURR_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]);ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;
    
    
    % Curr Corr
    IND = intersect(find(ISCORR_CURR(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_H);
    C_CURR_NUM=length(IND);
    subplot(2,2,3)
    C_CURR_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim; ylabel('Curr','fontsize',12);
    subplot(2,2,4)
    C_CURR_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;
    
    
    
    
    %MAIN TASK
    subplot(2,2,1)
    IND = intersect(CHANGE-15:CHANGE,IND_H);
    MT_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Prev','fontsize',12);
    subplot(2,2,2)
    MT_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim;
    
    subplot(2,2,3)
    IND = intersect(CHANGE-15:CHANGE,IND_H);
    MT_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Curr','fontsize',12);
    subplot(2,2,4)
    MT_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim;
    
    
    
    
    for i=1:4
        subplot(2,2,i)
        ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
        hold on;
        plot([0 0],ylim,'-k','linewidth',1);
    end
    
    
    suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12),'-',HAND));
    
    cd(Results_dir)
    filename = strcat('PSTH_Simple spike_DOSyDOS_TODO_',HAND);
    print(FF, '-dpdf', filename, '-r400')
   
    
    
    if ind_hand==1
        LEFT.WN_T = W_PREV_T;
        LEFT.CN_T = C_PREV_T;
        LEFT.WN_M = W_PREV_M;
        LEFT.CN_M = C_PREV_M;
        LEFT.NW_T = W_CURR_T;
        LEFT.NC_T = C_CURR_T;
        LEFT.NW_M = W_CURR_M;
        LEFT.NC_M = C_CURR_M;
    end
    
    if ind_hand==2
        RIGHT.WN_T = W_PREV_T;
        RIGHT.CN_T = C_PREV_T;
        RIGHT.WN_M = W_PREV_M;
        RIGHT.CN_M = C_PREV_M;
        RIGHT.NW_T = W_CURR_T;
        RIGHT.NC_T = C_CURR_T;
        RIGHT.NW_M = W_CURR_M;
        RIGHT.NC_M = C_CURR_M;
    end
    
    
    
end





%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






for index=1:4
    CW_DELTA.LEFT.statsWN(index) = NaN;
    CW_DELTA.LEFT.statsCN(index) = NaN;
    CW_DELTA.LEFT.statsNW(index) = NaN;
    CW_DELTA.LEFT.statsNC(index) = NaN;
end

for index=1:4
    CW_DELTA.RIGHT.statsWN(index) = NaN;
    CW_DELTA.RIGHT.statsCN(index) = NaN;
    CW_DELTA.RIGHT.statsNW(index) = NaN;
    CW_DELTA.RIGHT.statsNC(index) = NaN;
end


% ----------------------------------------------------------------------- %


for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        CW_DELTA.LEFT.statsWN(index) = nanmean(LEFT.WN_T(1,INDS));
        CW_DELTA.LEFT.statsCN(index) = nanmean(LEFT.CN_T(1,INDS));
        CW_DELTA.LEFT.statsNW(index) = NaN;
        CW_DELTA.LEFT.statsNC(index) = NaN;
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.LEFT.statsWN(index) = nanmean(LEFT.WN_T(1,INDS));
    CW_DELTA.LEFT.statsCN(index) = nanmean(LEFT.CN_T(1,INDS));
    CW_DELTA.LEFT.statsNW(index) = NaN;
    CW_DELTA.LEFT.statsNC(index) = NaN;
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.LEFT.statsWN(index) = NaN;
    CW_DELTA.LEFT.statsCN(index) = NaN;
    CW_DELTA.LEFT.statsNW(index) = nanmean(LEFT.NW_M(1,INDS));
    CW_DELTA.LEFT.statsNC(index) = nanmean(LEFT.NC_M(1,INDS));
end


% ----------------------------------------------------------------------- %


for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        CW_DELTA.RIGHT.statsWN(index) = nanmean(RIGHT.WN_T(1,INDS));
        CW_DELTA.RIGHT.statsCN(index) = nanmean(RIGHT.CN_T(1,INDS));
        CW_DELTA.RIGHT.statsNW(index) = NaN;
        CW_DELTA.RIGHT.statsNC(index) = NaN;
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.RIGHT.statsWN(index) = nanmean(RIGHT.WN_T(1,INDS));
    CW_DELTA.RIGHT.statsCN(index) = nanmean(RIGHT.CN_T(1,INDS));
    CW_DELTA.RIGHT.statsNW(index) = NaN;
    CW_DELTA.RIGHT.statsNC(index) = NaN;
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.RIGHT.statsWN(index) = NaN;
    CW_DELTA.RIGHT.statsCN(index) = NaN;
    CW_DELTA.RIGHT.statsNW(index) = nanmean(RIGHT.NW_M(1,INDS));
    CW_DELTA.RIGHT.statsNC(index) = nanmean(RIGHT.NC_M(1,INDS));
end





save(POP_file,'CW_DELTA','-append');
save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');




