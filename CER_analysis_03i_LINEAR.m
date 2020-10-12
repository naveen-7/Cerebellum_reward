

% function CER_analysis_03i_LINEAR
% written by naveen at cumc on 8/21/17


Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;


Type8 = NaN(length(Infos),1);
for i=2:length(Infos)
    if (Infos(i-1,10)==1)
        Type8(i,1)=1;
    end
end


Type3 = NaN(length(Infos),1);
for i=2:length(Infos)
    if (Infos(i-1,10)==2)
        Type3(i,1)=1;
    end
end




F = figure();


for index=1:4
    if DELTA.EpochFlag(index)==1
        
        index
        
        if index==1 EPOCH_TYPE = 'B'; end
        if index==2 EPOCH_TYPE = 'T'; end
        if index==3 EPOCH_TYPE = 'M'; end
        if index==4 EPOCH_TYPE = 'P'; end
        
        
        EPOCH_START = DELTA.START(index);
        EPOCH_END = DELTA.END(index);
        
        
        
        
        % FIRST THINGS ------------
        
        % Start = EPOCH_START;
        % End = EPOCH_END;
        
        if strcmp(EPOCH_TYPE,'P') | strcmp(EPOCH_TYPE,'M')
            Align_code = 11;
            ALIGN = 'M';
            Start = Start_M; End = End_M;
        end
        
        if strcmp(EPOCH_TYPE,'B') | strcmp(EPOCH_TYPE,'T')
            Align_code = 4;
            ALIGN = 'T';
            Start = Start_T; End = End_T;
        end
        
        
        TIME = Start:End;
        
        
        clear PSTH1 Matrix1 PSTH2 Matrix2 MAT
        for i = 1:size(Infos,1)  % Running for all trials in each case
            Sig = Spikes.S{i,1};
            PSTH1 = PSTH_ONE_n(Sig,Infos(i,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
            Matrix1(i,:) = nanmean(PSTH1(find(TIME==EPOCH_START+75):find(TIME==EPOCH_END-75)));
        end
        
        
        
        if index<4
            
            CORR_TRIALS1 = find(Type8(1:CHANGE)==1);
            CORR_VAL1 = Matrix1(CORR_TRIALS1,1);
            
            CORR_TRIALS2 = find(Type8(CHANGE+1:size(Infos,1))==1)+CHANGE+1-1;
            CORR_VAL2 = Matrix1(CORR_TRIALS2,1);
            
            CORR_TRIALS = [CORR_TRIALS1; CORR_TRIALS2];
            CORR_VAL = [CORR_VAL1; CORR_VAL2];
            
            
            WRNG_TRIALS1 = find(Type3(1:CHANGE)==1);
            WRNG_VAL1 = Matrix1(WRNG_TRIALS1,1);
            
            WRNG_TRIALS2 = find(Type3(CHANGE+1:size(Infos,1))==1)+CHANGE+1-1;
            WRNG_VAL2 = Matrix1(WRNG_TRIALS2,1);
            
            WRNG_TRIALS = [CORR_TRIALS1; CORR_TRIALS2];
            WRNG_VAL = [WRNG_VAL1; WRNG_VAL2];
            
            
            
            CHANGE_CORRECTED_CORR = CHANGE-length(find(Type3(1:CHANGE)==1));
            CHANGE_CORRECTED_WRNG = CHANGE-length(find(Type8(1:CHANGE)==1));
            
        end
        
        
        
        if index==4
            
            CORR_TRIALS1 = find(Infos(1:CHANGE,10)==1);
            CORR_VAL1 = Matrix1(CORR_TRIALS1,1);
            
            CORR_TRIALS2 = find(Infos(CHANGE+1:size(Infos,1),10)==1)+CHANGE+1-1;
            CORR_VAL2 = Matrix1(CORR_TRIALS2,1);
            
            CORR_TRIALS = [CORR_TRIALS1; CORR_TRIALS2];
            CORR_VAL = [CORR_VAL1; CORR_VAL2];
            
            WRNG_TRIALS1 = find(Infos(1:CHANGE,10)==2);
            WRNG_VAL1 = Matrix1(WRNG_TRIALS1,1);
            
            WRNG_TRIALS2 = find(Infos(CHANGE+1:size(Infos,1),10)==2)+CHANGE+1-1;
            WRNG_VAL2 = Matrix1(WRNG_TRIALS2,1);
            
            WRNG_TRIALS = [WRNG_TRIALS1; WRNG_TRIALS2];
            WRNG_VAL = [WRNG_VAL1; WRNG_VAL2];
            
            
            CHANGE_CORRECTED_CORR = CHANGE-length(find(Infos(1:CHANGE,10)==2));
            CHANGE_CORRECTED_WRNG = CHANGE-length(find(Infos(1:CHANGE,10)==1));
            
        end
             
        
        
        %           CV = smooth(CORR_VAL,0.3,'rloess');
        CV1_CORR = smooth(CORR_VAL1,0.3,'moving');
        CV2_CORR = smooth(CORR_VAL2,0.3,'moving');

        CV1_WRNG = smooth(WRNG_VAL1,0.3,'moving');
        CV2_WRNG = smooth(WRNG_VAL2,0.3,'moving');
        
        subplot(4,4,4+index)
        hold on;
        plot(1:length(CV1_CORR),CV1_CORR,'-B')
        plot(length(CV1_CORR)+1:length(CV2_CORR)+length(CV1_CORR),CV2_CORR,'-B')
        plot([CHANGE_CORRECTED_CORR CHANGE_CORRECTED_CORR],ylim,'-k')
        
        xlim([5 length(CV2_CORR)+length(CV1_CORR)-5]);
        set(gca,'linewidth',0.6,'fontsize',8,'fontweight','normal')
        
        
        
        subplot(4,4,8+index)
        hold on;
        plot(1:length(CV1_WRNG),CV1_WRNG,'-R')
        plot(length(CV1_WRNG)+1:length(CV2_WRNG)+length(CV1_WRNG),CV2_WRNG,'-R')
        plot([CHANGE_CORRECTED_WRNG CHANGE_CORRECTED_WRNG],ylim,'-k')
        
        xlim([1 length(CV2_WRNG)+length(CV1_WRNG)]);
        set(gca,'linewidth',0.6,'fontsize',8,'fontweight','normal')
        
        
        
        %     W_bin = 0.1; S_bin = 0.05;
        %     [x_trial, SIG_trial, per_corr] = LEARNING_CURVE_n(CORR_VAL, CHANGE_CORRECTED, Infos(CORR_TRIALS,:), W_bin, S_bin, 'p');
        %     figure();
        %     plot(x_trial,SIG_trial(:,1),'O')
        %     hold on;
        %     plot([CHANGE_CORRECTED CHANGE_CORRECTED],ylim,'-k')
        
        
        
        
        %         title(strcat(NOME(1:8),'-',NOME(10:12),'-CORR-STATE-',EPOCH_TYPE));
        %         cd(Results_dir)
        %         filename = strcat(NOME,'_CORR_STATE_Trial_by_Trail_',EPOCH_TYPE)
        %         print(F, '-dpdf', filename, '-r600')
        

        
    end
end






%% STATE CHANGE ANALYSIS PROPER





% xlabel('trial number')
% ylabel('sp/s')







suptitle(strcat('STATE-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = strcat(NOME,'_STATE_TRIAL_BY_TRIAL_CW');
print(F, '-dpdf', filename, '-r400')









%%%% reaction time analysis -------------------------------

F = figure();

IND_CORR = find(Infos(CHANGE:CHANGE+19,10)==1)+CHANGE-1;
IND_WRNG = find(Infos(CHANGE:CHANGE+19,10)==2)+CHANGE-1;

clear RT

RT.CORR_SMALL = NaN(1,4);
RT.CORR_LARGE = NaN(1,4);

RT.WRNG_SMALL = NaN(1,4);
RT.WRNG_LARGE = NaN(1,4);
        

for index=1:4
    if DELTA.EpochFlag(index)==1
        
        index
        
        if index==1 EPOCH_TYPE = 'B'; end
        if index==2 EPOCH_TYPE = 'T'; end
        if index==3 EPOCH_TYPE = 'M'; end
        if index==4 EPOCH_TYPE = 'P'; end
        
        
        EPOCH_START = DELTA.START(index);
        EPOCH_END = DELTA.END(index);
        
        
        if strcmp(EPOCH_TYPE,'P') | strcmp(EPOCH_TYPE,'M')
            Align_code = 11;
            ALIGN = 'M';
            Start = Start_M; End = End_M;
        end
        
        if strcmp(EPOCH_TYPE,'B') | strcmp(EPOCH_TYPE,'T')
            Align_code = 4;
            ALIGN = 'T';
            Start = Start_T; End = End_T;
        end
        
        
        TIME = Start:End;
        
        
        clear PSTH1 Matrix1 PSTH2 Matrix2 MAT
        for i = 1:size(Infos,1)  % Running for all trials in each case
            Sig = Spikes.S{i,1};
            PSTH1 = PSTH_ONE_n(Sig,Infos(i,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
            Matrix1(i,:) = nanmean(PSTH1(find(TIME==EPOCH_START+75):find(TIME==EPOCH_END-75)));
        end
        
        
        ALL_CORR = Matrix1(IND_CORR);
        ALL_WRNG = Matrix1(IND_WRNG);
        
        RT_CORR = Infos(IND_CORR,14);
        RT_WRNG = Infos(IND_WRNG,14);
        
        
        % CLASSIFY RT
        IND_CORR_SMALL = find(RT_CORR<= quantile(RT_CORR,0.5));
        IND_CORR_LARGE = find(RT_CORR>= quantile(RT_CORR,0.5));
        
        IND_WRNG_SMALL = find(RT_WRNG<= quantile(RT_WRNG,0.5));
        IND_WRNG_LARGE = find(RT_WRNG>= quantile(RT_WRNG,0.5));
        
        
        RT.CORR_SMALL(index) = nanmean(ALL_CORR(IND_CORR_SMALL));
        RT.CORR_LARGE(index) = nanmean(ALL_CORR(IND_CORR_LARGE));
        
        RT.WRNG_SMALL(index) = nanmean(ALL_WRNG(IND_WRNG_SMALL));
        RT.WRNG_LARGE(index) = nanmean(ALL_WRNG(IND_WRNG_LARGE));
        
        
        
        subplot(4,4,index*4-3)
        hold on;
        bar([1 2],[RT.CORR_SMALL(index) RT.CORR_LARGE(index)],'b')
        bar([4 5],[RT.WRNG_SMALL(index) RT.WRNG_LARGE(index)],'r')
        xlim([0 6])
        box off;
        
    end
end



suptitle(strcat('RT EFFECT-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = strcat(NOME,'_RT_EFFECT');
print(F, '-dpdf', filename, '-r400')



save(POP_file,'RT','-append');
save(MERGE_file,'RT','-append');
save(ALLCELLS_file,'RT','-append');





% end