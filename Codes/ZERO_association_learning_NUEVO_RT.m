
%% function ZERO_association_learning_NUEVO_RT

% Written by naveen on 2/15/18 at cumc

% Take 5 trials before change and 5 trials after change
% take 5 random trials in bef and 5 random trials in aft and take mean
% compute distance for the mean
% repeat 250 times
% calculate the 'maximum allowed variance' in the original data in the same way




RT = Infos(:,14);
RT1 = Infos(:,14);
RT = RT(CHANGE:CHANGE+30);
RTcomp{1,:} = (find(RT<median(RT)))+CHANGE-1;
RTcomp{2,:} = (find(RT>median(RT)))+CHANGE-1;

clear RTVAL;
RTVAL(:,1) = RT1(CHANGE-15:CHANGE-1);
RTVAL(:,2) = RT(find(RT<median(RT)));
RTVAL(:,3) = RT(find(RT>median(RT)));

Ass_Learn.RTVAL = RTVAL;



% % % 
% % % 
% % % 
% % % 
% % % 
% % % NUM = 250;
% % % 
% % % Start_T = -450;
% % % End_T  =1050;
% % % 
% % % StartLIM_T = -400; EndLIM_T = 400;
% % % TIME_T = Start_T:End_T;
% % % 
% % % Start_M = -750;
% % % End_M  =750;
% % % StartLIM_M = -200; EndLIM_M = 600;
% % % TIME_M = Start_M:End_M;
% % % 
% % % 
% % % num_trials = 15;
% % % RANGE = CHANGE-num_trials:CHANGE+num_trials-1;
% % % Start_time = Start_M;
% % % End_time = End_M;
% % % Sigma = 30;
% % % 
% % % GROUP1 = CHANGE-num_trials:CHANGE-1;
% % % 
% % % 
% % % 
% % % 
% % % for kkk=1:NUM
% % % 
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%% within condition %%%%%%%%%%%%%%
% % % for i=1:2
% % %     IND = randperm(15);
% % %     IND = IND(1:5);
% % %     
% % %     % part 1----------------------------------------
% % %     clear trials_corr PSTH_corr trials_wrng P1 ind TIME
% % %     Start_time = -600; End_time = 400; TIME = Start_time:End_time;
% % %     START = -400; END = 200;
% % %     P1 = PSTH_RETURN_n(Spikes.S(GROUP1(IND)),Infos(GROUP1(IND),4),Start_time,End_time,Sigma);
% % %     
% % %     % part 2----------------------------------------
% % %     clear trials_corr PSTH_corr trials_wrng P2 ind TIME
% % %     Start_time = -400; End_time = 800; TIME = Start_time:End_time;
% % %     START = -200; END = 600;
% % %     P2 = PSTH_RETURN_n(Spikes.S(GROUP1(IND)),Infos(GROUP1(IND),11),Start_time,End_time,Sigma);
% % %     
% % %     P_within = [P1 P2];
% % %     P_within_avg(i,:) = nanmean(P_within);
% % % end
% % % 
% % % WITHIN(kkk) =  sqrt(nansum((P_within_avg(1,:)-P_within_avg(2,:)).^2));
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%% across condition %%%%%%%%%%%%%%
% % % RT_USE = RTcomp{1,:};
% % % GROUP2 = RT_USE;
% % % %%%%%% step 1: OT CONDITION
% % % IND = randperm(15);
% % % IND = IND(1:5);
% % % 
% % % % part 1----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P1 ind TIME
% % % Start_time = -600; End_time = 400; TIME = Start_time:End_time;
% % % START = -400; END = 200;
% % % P1 = PSTH_RETURN_n(Spikes.S(GROUP1(IND)),Infos(GROUP1(IND),4),Start_time,End_time,Sigma);
% % % 
% % % % part 2----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P2 ind TIME
% % % Start_time = -400; End_time = 800; TIME = Start_time:End_time;
% % % START = -200; END = 600;
% % % P2 = PSTH_RETURN_n(Spikes.S(GROUP1(IND)),Infos(GROUP1(IND),11),Start_time,End_time,Sigma);
% % % 
% % % P_OT = nanmean([P1 P2]);
% % % 
% % % 
% % % %%%%%% step 2: N CONDITION
% % % IND = randperm(15);
% % % IND = IND(1:5);
% % % 
% % % % part 1----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P1 ind TIME
% % % Start_time = -600; End_time = 400; TIME = Start_time:End_time;
% % % START = -400; END = 200;
% % % P1 = PSTH_RETURN_n(Spikes.S(GROUP2(IND)),Infos(GROUP2(IND),4),Start_time,End_time,Sigma);
% % % 
% % % % part 2----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P2 ind TIME
% % % Start_time = -400; End_time = 800; TIME = Start_time:End_time;
% % % START = -200; END = 600;
% % % P2 = PSTH_RETURN_n(Spikes.S(GROUP2(IND)),Infos(GROUP2(IND),11),Start_time,End_time,Sigma);
% % % 
% % % P_N = nanmean([P1 P2]);
% % % 
% % % 
% % % ACROSS_FAST(kkk) =  sqrt(nansum((P_OT-P_N).^2));
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%% across condition %%%%%%%%%%%%%%
% % % RT_USE = RTcomp{2,:};
% % % GROUP2 = RT_USE;
% % % %%%%%% step 1: OT CONDITION
% % % IND = randperm(15);
% % % IND = IND(1:5);
% % % 
% % % % part 1----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P1 ind TIME
% % % Start_time = -600; End_time = 400; TIME = Start_time:End_time;
% % % START = -400; END = 200;
% % % P1 = PSTH_RETURN_n(Spikes.S(GROUP1(IND)),Infos(GROUP1(IND),4),Start_time,End_time,Sigma);
% % % 
% % % % part 2----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P2 ind TIME
% % % Start_time = -400; End_time = 800; TIME = Start_time:End_time;
% % % START = -200; END = 600;
% % % P2 = PSTH_RETURN_n(Spikes.S(GROUP1(IND)),Infos(GROUP1(IND),11),Start_time,End_time,Sigma);
% % % 
% % % P_OT = nanmean([P1 P2]);
% % % 
% % % 
% % % %%%%%% step 2: N CONDITION
% % % IND = randperm(15);
% % % IND = IND(1:5);
% % % 
% % % % part 1----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P1 ind TIME
% % % Start_time = -600; End_time = 400; TIME = Start_time:End_time;
% % % START = -400; END = 200;
% % % P1 = PSTH_RETURN_n(Spikes.S(GROUP2(IND)),Infos(GROUP2(IND),4),Start_time,End_time,Sigma);
% % % 
% % % % part 2----------------------------------------
% % % clear trials_corr PSTH_corr trials_wrng P2 ind TIME
% % % Start_time = -400; End_time = 800; TIME = Start_time:End_time;
% % % START = -200; END = 600;
% % % P2 = PSTH_RETURN_n(Spikes.S(GROUP2(IND)),Infos(GROUP2(IND),11),Start_time,End_time,Sigma);
% % % 
% % % P_N = nanmean([P1 P2]);
% % % 
% % % 
% % % ACROSS_SLOW(kkk) =  sqrt(nansum((P_OT-P_N).^2));
% % % 
% % % 
% % % end
% % % 
% % % 
% % % ACROSS_color = [33 178 128]/255;
% % % WITHIN_color = [102 102 102]/255;
% % % 
% % % 
% % % BIN_W = 40;
% % % F = figure;
% % % subplot(2,2,1,'box','off','tickdir','out','Linewidth',0.3,'FontSize',7)
% % % hold on;
% % % h = histogram(WITHIN,'BinWidth',BIN_W);
% % % h.FaceColor = WITHIN_color;
% % % h = histogram(ACROSS_FAST,'BinWidth',BIN_W);
% % % h.FaceColor = ACROSS_color;
% % % box off;
% % % xlabel('rms distance(a.u.)');
% % % ylabel('frequency');
% % % xlim([0 1000]); ylim([0 100]);
% % % title('fast');
% % % 
% % % 
% % % subplot(2,2,3,'box','off','tickdir','out','Linewidth',0.3,'FontSize',7)
% % % hold on;
% % % h = histogram(WITHIN,'BinWidth',BIN_W);
% % % h.FaceColor = WITHIN_color;
% % % h = histogram(ACROSS_SLOW,'BinWidth',BIN_W);
% % % h.FaceColor = ACROSS_color;
% % % box off;
% % % xlabel('rms distance(a.u.)');
% % % ylabel('frequency');
% % % xlim([0 1000]); ylim([0 100]);
% % % title('slow');
% % % 
% % %  cd(Results_dir)
% % %         filename = strcat(NOME,'_Association_Learning_MRD_STATS_RT');
% % %         print(F, '-dpdf', filename, '-r400')
% % %     
% % % 
% % % 
% % % if ( lillietest(WITHIN)==0 & lillietest(ACROSS_FAST)==0 & lillietest(ACROSS_SLOW)==0 )
% % %     Ass_Learn.WITHIN = nanmean(WITHIN);
% % %     Ass_Learn.ACROSSF = nanmean(ACROSS_FAST);
% % %     Ass_Learn.ACROSSS = nanmean(ACROSS_SLOW);
% % % else
% % %     Ass_Learn.WITHIN = nanmedian(WITHIN);
% % %     Ass_Learn.ACROSSF = nanmedian(ACROSS_FAST);
% % %     Ass_Learn.ACROSSS = nanmedian(ACROSS_SLOW);
% % % end
% % % 
% % % 
% % % Ass_Learn.P_WITHIN_ACROSSF = ttest_NN(WITHIN,ACROSS_FAST);
% % % Ass_Learn.P_WITHIN_ACROSSS = ttest_NN(WITHIN,ACROSS_SLOW);



if strcmp(POP_file(1),'C') POP_file(1)='E'; end
if strcmp(MERGE_file(1),'C') MERGE_file(1)='E'; end
if strcmp(Results_dir(1),'C') Results_dir(1)='E'; end
if strcmp(ALLCELLS_file(1),'C') ALLCELLS_file(1)='E'; end

save(POP_file,'Ass_Learn','-append');
save(MERGE_file,'Ass_Learn','-append');
save(ALLCELLS_file,'Ass_Learn','-append');




% end