%% function CER_analysis_03m_LINEAR
% does sort of a bootstrap analysis to proove that delta epoch is real ------------------
% written by naveen at JLG on 10/17/18

cycles = 500;
NUM = 20;

%% %%%%%%% SHOESTRAPPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIG = 40;

Start_T = -650;
End_T  =1050;
StartLIM_T = -500; EndLIM_T = 900;

Start_M = -750;
End_M  =950;
StartLIM_M = -500; EndLIM_M = 900;

DELTA_PROOF.NULL = NaN(1,4);
DELTA_PROOF.SHUF = NaN(1,4);
DELTA_PROOF.p = NaN(1,4);

LOC = find(DELTA.EpochFlag==1);

% if length(loc)==1 %%%%%%%%%%%%%%%%% ONLY DO SINGLE DELTA neurons

for ii=1:(length(LOC))

    loc = LOC(ii)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:1 initialize everything according to delta location
    if loc<=2
        Align = 4;
        Start = Start_T; End = End_T; StartLIM = StartLIM_T; EndLIM = EndLIM_T;
        psthtemp=PSTH_RETURN_n(Spikes.S(CHANGE:CHANGE+NUM-1,1),Infos(CHANGE:CHANGE+NUM-1,11),Start,End,SIG);
        time = Start:End;
        SDF_tru = psthtemp(:,find(time==StartLIM): find(time==EndLIM));
        D_srt = DELTA.START(loc);
        D_end = DELTA.END(loc);
        time_tru = StartLIM:EndLIM;
    end
    if loc==3|loc==4
        Align = 11;
        Start = Start_M; End = End_M; StartLIM = StartLIM_M; EndLIM = EndLIM_M;
        psthtemp=PSTH_RETURN_n(Spikes.S(CHANGE:CHANGE+NUM-1,1),Infos(CHANGE:CHANGE+NUM-1,11),Start,End,SIG);
        time = Start:End;
        SDF_tru = psthtemp(:,find(time==StartLIM): find(time==EndLIM));
        D_srt = DELTA.START(loc);
        D_end = DELTA.END(loc);
        time_tru = StartLIM:EndLIM;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:2 create matrix of PSTH and trial outcomes amd time
    
    if loc<=3
        VAL_tru = Infos(CHANGE-1:CHANGE+NUM-2,10);
    end
    
    if loc==4
        VAL_tru = Infos(CHANGE:CHANGE+NUM-1,10);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:3 do shoestrap for null data
    
    for i=1:cycles
        SDF_use = SDF_tru;
        VAL_use = VAL_tru;
        SHUFF_IND = randperm(20);
        SDF_use = SDF_use(SHUFF_IND,:);
        VAL_use = VAL_use(SHUFF_IND,:);
        
        corr_ind = find(VAL_use==1,5);
        wrng_ind = find(VAL_use==2,5);
        
        PLOT1(1,:)=nanmean(SDF_use(corr_ind,:));
        PLOT1(2,:)=nanmean(SDF_use(wrng_ind,:));
        
        
        EPOCH = find(time_tru==D_srt):find(time_tru==D_end);
        clear COMPARE;
         COMPARE(1,:)=nanmean(SDF_use(corr_ind,EPOCH));
         COMPARE(2,:)=nanmean(SDF_use(wrng_ind,EPOCH));
       
        NULL(i,1) = sqrt(sum((COMPARE(1,:)-COMPARE(2,:)).^2));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:4 do shoestrap for shuffle data
    
    for i=1:cycles
        SDF_use = SDF_tru;
        VAL_use = VAL_tru;
        SHUFF_IND = randperm(20);
        %SDF_use = SDF_use(SHUFF_IND,:);
        VAL_use = VAL_use(SHUFF_IND,:);
        
        corr_ind = find(VAL_use==1,5);
        wrng_ind = find(VAL_use==2,5);
        
        PLOT2(1,:)=nanmean(SDF_use(corr_ind,:));
        PLOT2(2,:)=nanmean(SDF_use(wrng_ind,:));
        
        
        EPOCH = find(time_tru==D_srt):find(time_tru==D_end);
        clear COMPARE;
         COMPARE(1,:)=nanmean(SDF_use(corr_ind,EPOCH));
         COMPARE(2,:)=nanmean(SDF_use(wrng_ind,EPOCH));
         
        SHUFF(i,1) = sqrt(sum((COMPARE(1,:)-COMPARE(2,:)).^2));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:4 plotting nd stuff
    
    F = figure()
    subplot(2,2,1)
    hold on;
    plot(time_tru,nanmean(SDF_tru(find(VAL_tru==1),:)),'-b')
    plot(time_tru,nanmean(SDF_tru(find(VAL_tru==2),:)),'-r')
    
    subplot(2,2,2)
    hold on;
    plot(time_tru(EPOCH),nanmean(SDF_tru(find(VAL_tru==1),EPOCH)),'-b')
    plot(time_tru(EPOCH),nanmean(SDF_tru(find(VAL_tru==2),EPOCH)),'-r')
    
    
%     COM(1,:)=nanmean(SDF_tru(find(VAL_tru==1),EPOCH));
%     COM(2,:)=nanmean(SDF_tru(find(VAL_tru==2),EPOCH));
%     COMT = mat2gray([COM(1,:); COM(2,:)]);
%     COM(1,:) = COMT(1,:); COM(2,:) = COMT(2,:);
%     
%     subplot(2,2,3)
%     hold on;
%     plot(time_tru(EPOCH),(COM(1,:)),'-b')
%     plot(time_tru(EPOCH),(COM(2,:)),'-r')
    
    
    subplot(2,2,4)
    hold on;
    NULL_color = [209 173 21]/255;
    SHUFF_color = [102 102 102]/255;
    BIN_W = 40;
    h = histogram(NULL,'BinWidth',BIN_W);
    h.FaceColor = NULL_color;
    h = histogram(SHUFF,'BinWidth',BIN_W);
    h.FaceColor = SHUFF_color;
    box off;
    xlabel('rms distance(a.u.)');
    ylabel('frequency');
    YLIM = ylim;
    p = stats_test_n(NULL,SHUFF);
    plot([nanmean(SHUFF) nanmean(NULL)],[YLIM(2) YLIM(2)],'-k')
    text(nanmean(SHUFF)+0.25*(nanmean(NULL)-nanmean(SHUFF)), YLIM(2), star_n(p));
    
    
    cd(Results_dir)
    filename = strcat(NOME,'DELTA_MRD');
    print(F, '-dpdf', filename, '-r400')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:5 retrun these stuff
    DELTA_PROOF.NULL(LOC(ii)) = nanmean(NULL);
    DELTA_PROOF.SHUF(LOC(ii)) = nanmean(SHUFF);
    DELTA_PROOF.p(LOC(ii)) = stats_test_n(NULL,SHUFF);
    
    
    
    
    
end
    
    
    
%     
%   F =  figure()
%    subplot(2,2,1)
%    hold on;
%    plot(time_tru,PLOT1(1,:),'-b')
%    plot(time_tru,PLOT1(2,:),'-r')
%    ylim([0 80]);
%    plot([StartLIM StartLIM],[0 20],'-k')
%    plot([StartLIM StartLIM+200],[0 0],'-k')
%    plot([0 0],[0 0],'ok','markersize',3)
%    axis off;
%     
%    subplot(2,2,2)
%    hold on;
%    plot(time_tru,PLOT2(1,:),'-b')
%    plot(time_tru,PLOT2(2,:),'-r')
%    ylim([0 80]);
%    plot([StartLIM StartLIM],[0 20],'-k')
%    plot([StartLIM StartLIM+200],[0 0],'-k')
%    plot([0 0],[0 0],'ok','markersize',3)
%    axis off;
%     
%        
%     cd(Results_dir)
%     filename = strcat(NOME,'DELTA_MRD_example');
%     print(F, '-dpdf', filename, '-r400')
%      
   
    
    
%  end






save(POP_file,'DELTA_PROOF','-append');
save(MERGE_file,'DELTA_PROOF','-append');
save(ALLCELLS_file,'DELTA_PROOF','-append');
