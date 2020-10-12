
%% function CER_analysis_A05_ADVANCED_PLUGIN6
% doing the memory analysis for each symbol ------------------

% written by naveen at cumc on 8/12/17
%

load(POP_file)

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




%%%%%%%%%%%%%%%%%%%%%%% stimuli 

if isnan(CHANGE)
    CHANGE=1;
end
  
temp1 = NaN(size(Infos,1),1);
Stim  = NaN(size(Infos,1),1);
flag  = zeros(4,1);
temp2 = 0;
temp3 = NaN(size(Infos,1),1);
temps = NaN(4,1);

for i=CHANGE:CHANGE+20

    if Infos(i,5)==160
        temp1(i) = (4*(Infos(i,6))) + (9*(Infos(i,7))) + (16*(Infos(i,8)));
    end

    if Infos(i,5)==160
        if (temp1(1) == temp1(i))
            Stim(i)=1; else Stim(i)=2; end
    end


    if (Infos(i,5)>160) && (flag(1)==0)
        flag(1)=1;
        temps(1) = Infos(i,5);
    end

    if (Infos(i,5)>160) && (flag(1)==1)
        temp3(i) = Infos(i,5);

        if (temp3(i) ~= temps(1)) && (flag(2)==0)
            temps(2)=temp3(i);
            flag(2)=1;
        end

        if (temp3(i) ~= temps(1)) && (temp3(i) ~= temps(2)) && (flag(2)==1) && (flag(3)==0)
            temps(3)=temp3(i);
            flag(3)=1;
        end

        if (temp3(i) ~= temps(1)) && (temp3(i) ~= temps(2)) && (temp3(i) ~= temps(3)) && (flag(3)==1) && (flag(4)==0)
            temps(4)=temp3(i);
            flag(4)=1;
        end

    end


    if (Infos(i,5)>160) && (flag(1)==1)
        if (temp3(i) == temps(1))
            Stim(i)=3;
        elseif (temp3(i) == temps(2))
            Stim(i)=4;
        elseif (temp3(i) == temps(3))
            Stim(i)=5;
        elseif (temp3(i) == temps(4))
            Stim(i)=6;
        end
    end
    

end  
    


if nanmean(~(Infos(:,5)==160))==1
    SYM = unique(Infos(CHANGE:CHANGE+10,5));
    IND_ONE = find(Infos(:,5)==SYM(1));
    IND_TWO = find(Infos(:,5)==SYM(2));
end













for ind_targ=1:2
    
    if  ind_targ==1
        IND_T=IND_ONE;
        HAND = 'ONE';
    end
    
    if  ind_targ==2
        IND_T=IND_TWO;
        HAND = 'TWO';
    end
    
    
    FF = figure;
    SIG =40;
    
    
    % Prev Wrng
    IND = intersect(find(ISWRNG_PREV(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_T);
    W_PREV_NUM=length(IND);
    subplot(2,2,1)
    W_PREV_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
    subplot(2,2,2)
    W_PREV_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim;
    
    
    % Prev Corr
    IND = intersect(find(ISCORR_PREV(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_T);
    C_PREV_NUM=length(IND);
    subplot(2,2,1)
    C_PREV_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim; ylabel('Prev','fontsize',12);
    subplot(2,2,2)
    C_PREV_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;
    
    
    % Curr Wrng
    IND = intersect(find(ISWRNG_CURR(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_T);
    W_CURR_NUM=length(IND);
    subplot(2,2,3)
    W_CURR_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
    subplot(2,2,4)
    W_CURR_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]);ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;
    
    
    % Curr Corr
    IND = intersect(find(ISCORR_CURR(CHANGE:CHANGE+19,1)==1)+CHANGE-1,IND_T);
    C_CURR_NUM=length(IND);
    subplot(2,2,3)
    C_CURR_T = PSTH_n(SIGNAL(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim; ylabel('Curr','fontsize',12);
    subplot(2,2,4)
    C_CURR_M = PSTH_n(SIGNAL(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
    set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;
    
    
    
    suptitle(strcat('PSTH-SS-',NOME));
    
    cd(Results_dir)
    filename = strcat('PSTH_Simple spike_DOSyDOS_TODO_',HAND);
    print(FF, '-dpdf', filename, '-r400')
   
    
    
    if ind_targ==1
        ONE.WN_T = W_PREV_T;
        ONE.CN_T = C_PREV_T;
        ONE.WN_M = W_PREV_M;
        ONE.CN_M = C_PREV_M;
        ONE.NW_T = W_CURR_T;
        ONE.NC_T = C_CURR_T;
        ONE.NW_M = W_CURR_M;
        ONE.NC_M = C_CURR_M;
    end
    
    if ind_targ==2
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

if CHANGE==1
    CHANGE=NaN;
end



