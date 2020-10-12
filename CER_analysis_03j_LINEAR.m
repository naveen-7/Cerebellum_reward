
%% function CER_analysis_03j_LINEAR
% doing the memory analysis for each symbol ------------------

% written by naveen at cumc on 9/4/17
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



SYM = unique(Infos(CHANGE:CHANGE+10,5));

% SYMBOLS seperated -------------

IND_ONE = find(Infos(:,5)==SYM(1));
IND_TWO = find(Infos(:,5)==SYM(2));

for ind_symb=1:2
    
    if  ind_symb==1
        IND_S=IND_ONE;
        SYMBOL = 'ONE';
    end
    
    if  ind_symb==2
        IND_S=IND_TWO;
        SYMBOL = 'TWO';
    end
    
    
    FF = figure;
    SIG =40;
    
    
    % Prev Wrng
    IND = intersect(find(ISWRNG_PREV(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_S);
    W_PREV_NUM=length(IND);
    subplot(2,2,1)
    W_PREV_T = PSTHe_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
    subplot(2,2,2)
    W_PREV_M = PSTHe_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim;
    
    
    % Prev Corr
    IND = intersect(find(ISCORR_PREV(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_S);
    C_PREV_NUM=length(IND);
    subplot(2,2,1)
    C_PREV_T = PSTHe_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim; ylabel('Prev','fontsize',12);
    subplot(2,2,2)
    C_PREV_M = PSTHe_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;
    
    
    % Curr Wrng
    IND = intersect(find(ISWRNG_CURR(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_S);
    W_CURR_NUM=length(IND);
    subplot(2,2,3)
    W_CURR_T = PSTHe_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
    subplot(2,2,4)
    W_CURR_M = PSTHe_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]);ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;
    
    
    % Curr Corr
    IND = intersect(find(ISCORR_CURR(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_S);
    C_CURR_NUM=length(IND);
    subplot(2,2,3)
    C_CURR_T = PSTHe_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim; ylabel('Curr','fontsize',12);
    subplot(2,2,4)
    C_CURR_M = PSTHe_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;
    
    
    
    
    %MAIN TASK
    subplot(2,2,1)
    IND = intersect(CHANGE-15:CHANGE,IND_S);
    MT_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Prev','fontsize',12);
    subplot(2,2,2)
    MT_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim;
    
    subplot(2,2,3)
    IND = intersect(CHANGE-15:CHANGE,IND_S);
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
    
    
    suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12),'-',SYMBOL));
    
    cd(Results_dir)
    filename = strcat('PSTH_Simple spike_DOSyDOS_TODO_',SYMBOL);
    print(FF, '-dpdf', filename, '-r400')
   
    
    
    if ind_symb==1
        ONE.WN_T = W_PREV_T;
        ONE.CN_T = C_PREV_T;
        ONE.WN_M = W_PREV_M;
        ONE.CN_M = C_PREV_M;
        ONE.NW_T = W_CURR_T;
        ONE.NC_T = C_CURR_T;
        ONE.NW_M = W_CURR_M;
        ONE.NC_M = C_CURR_M;
    end
    
    if ind_symb==2
        TWO.WN_T = W_PREV_T;
        TWO.CN_T = C_PREV_T;
        TWO.WN_M = W_PREV_M;
        TWO.CN_M = C_PREV_M;
        TWO.NW_T = W_CURR_T;
        TWO.NC_T = C_CURR_T;
        TWO.NW_M = W_CURR_M;
        TWO.NC_M = C_CURR_M;
    end
    
    
    
end





%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






for index=1:4
    CW_DELTA.ONE.statsWN(index) = NaN;
    CW_DELTA.ONE.statsCN(index) = NaN;
    CW_DELTA.ONE.statsNW(index) = NaN;
    CW_DELTA.ONE.statsNC(index) = NaN;
end

for index=1:4
    CW_DELTA.TWO.statsWN(index) = NaN;
    CW_DELTA.TWO.statsCN(index) = NaN;
    CW_DELTA.TWO.statsNW(index) = NaN;
    CW_DELTA.TWO.statsNC(index) = NaN;
end


% ----------------------------------------------------------------------- %


for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        CW_DELTA.ONE.statsWN(index) = nanmean(ONE.WN_T(1,INDS));
        CW_DELTA.ONE.statsCN(index) = nanmean(ONE.CN_T(1,INDS));
        CW_DELTA.ONE.statsNW(index) = NaN;
        CW_DELTA.ONE.statsNC(index) = NaN;
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.ONE.statsWN(index) = nanmean(ONE.WN_T(1,INDS));
    CW_DELTA.ONE.statsCN(index) = nanmean(ONE.CN_T(1,INDS));
    CW_DELTA.ONE.statsNW(index) = NaN;
    CW_DELTA.ONE.statsNC(index) = NaN;
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.ONE.statsWN(index) = NaN;
    CW_DELTA.ONE.statsCN(index) = NaN;
    CW_DELTA.ONE.statsNW(index) = nanmean(ONE.NW_M(1,INDS));
    CW_DELTA.ONE.statsNC(index) = nanmean(ONE.NC_M(1,INDS));
end


% ----------------------------------------------------------------------- %


for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
        CW_DELTA.TWO.statsWN(index) = nanmean(TWO.WN_T(1,INDS));
        CW_DELTA.TWO.statsCN(index) = nanmean(TWO.CN_T(1,INDS));
        CW_DELTA.TWO.statsNW(index) = NaN;
        CW_DELTA.TWO.statsNC(index) = NaN;
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.TWO.statsWN(index) = nanmean(TWO.WN_T(1,INDS));
    CW_DELTA.TWO.statsCN(index) = nanmean(TWO.CN_T(1,INDS));
    CW_DELTA.TWO.statsNW(index) = NaN;
    CW_DELTA.TWO.statsNC(index) = NaN;
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
    CW_DELTA.TWO.statsWN(index) = NaN;
    CW_DELTA.TWO.statsCN(index) = NaN;
    CW_DELTA.TWO.statsNW(index) = nanmean(TWO.NW_M(1,INDS));
    CW_DELTA.TWO.statsNC(index) = nanmean(TWO.NC_M(1,INDS));
end





save(POP_file,'CW_DELTA','-append');
save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');




