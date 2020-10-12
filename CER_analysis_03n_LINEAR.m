
%   CER_analysis_03n_LINEAR %%%% ONLY delta epoch provides a cognitive error


del_epo = find(DELTA.EpochFlag==1);
del_num = length(find(DELTA.EpochFlag==1));
Sigma = 40;

MEAN = NaN(4,4); MED = NaN(4,4); STD = NaN(4,4); ACTUAL = NaN(4,4);

Start_time_T = -500;  End_time_T = 900; TIME_T = Start_time_T:End_time_T;
Start_time_M = -500;  End_time_M = 900; TIME_M = Start_time_M:End_time_M;



for i=1:del_num
    
    F = figure()
    
    for ite = 1:4
        
        if ite==1 TRIALS = 2:CHANGE-1;           end
        if ite==2 TRIALS = CHANGE:CHANGE+20;     end
        if ite==3 TRIALS = CHANGE+21:LEARNT-10;  end
        if ite==4 TRIALS = LEARNT-11:size(Infos,1)-1; end
        
        
        
        DEL_epoch = del_epo(i);
        clear Corr_trials Wrng_trials
        
        if DEL_epoch == 4
            %         Align_code = 11;
            Corr_trials = find(Infos(TRIALS,10)==1)+TRIALS(1)-1;
            Wrng_trials = find(Infos(TRIALS,10)==2)+TRIALS(1)-1;
        end
        if DEL_epoch == 3
            %         Align_code = 11;
            Corr_trials = find(Infos(TRIALS,10)==1)+TRIALS(1);
            Wrng_trials = find(Infos(TRIALS,10)==2)+TRIALS(1);
        end
        if DEL_epoch == 1 | DEL_epoch == 2
            %         Align_code = 4;
            Corr_trials = find(Infos(TRIALS,10)==1)+TRIALS(1);
            Wrng_trials = find(Infos(TRIALS,10)==2)+TRIALS(1);
        end
        
        
        
        subplot(4,4,4*(ite-1)+1)
        PSTH_corr_T = PSTH_n(Spikes.S(Corr_trials),Infos(Corr_trials,4),Start_time_T,End_time_T,Sigma,[0 0 1],1,0);
        PSTH_wrng_T = PSTH_n(Spikes.S(Wrng_trials),Infos(Wrng_trials,4),Start_time_T,End_time_T,Sigma,[1 0 0],1,0);
        
        subplot(4,4,4*(ite-1)+2)
        PSTH_corr_M = PSTH_n(Spikes.S(Corr_trials),Infos(Corr_trials,11),Start_time_M,End_time_M,Sigma,[0 0 1],1,0);
        PSTH_wrng_M = PSTH_n(Spikes.S(Wrng_trials),Infos(Wrng_trials,11),Start_time_M,End_time_M,Sigma,[1 0 0],1,0);
        
        
        %%%%% cortando-----------------
        TIME_t = -400:400; TIME_m = -400:800;
        
        PSTH_corr_T = PSTH_corr_T(1,find(TIME_t(1)<=TIME_T & TIME_T<=TIME_t(end)));
        PSTH_wrng_T = PSTH_wrng_T(1,find(TIME_t(1)<=TIME_T & TIME_T<=TIME_t(end)));
        PSTH_corr_M = PSTH_corr_M(1,find(TIME_m(1)<=TIME_T & TIME_T<=TIME_m(end)));
        PSTH_wrng_M = PSTH_wrng_M(1,find(TIME_m(1)<=TIME_T & TIME_T<=TIME_m(end)));
        
        
        %%%%% GENERATE RANDOM LIST OF 200 ms time inervals
        
        for ii=1:200
            if rand(1)<=0.5
                time = TIME_t;
                EPOCHS(ii,3) = 4;
            else
                time = TIME_m;
                EPOCHS(ii,3) = 11;
            end
            EPOCHS(ii,1) = time(randi(length(time)-200));
            EPOCHS(ii,2) = EPOCHS(ii,1)+200;
        end
        
        
        %%%%% find the activity in the generated time bins
        
        for ii=1:200
            if  EPOCHS(ii,3) == 4
                DIFF(ii,1) = abs( nanmean(PSTH_corr_T(find( EPOCHS(ii,1)<=TIME_t & TIME_t<=EPOCHS(ii,2) ))) ...
                    - nanmean(PSTH_wrng_T(find( EPOCHS(ii,1)<=TIME_t & TIME_t<=EPOCHS(ii,2) ))) );
            end
            if  EPOCHS(ii,3) == 11
                DIFF(ii,1) = abs( nanmean(PSTH_corr_M(find( EPOCHS(ii,1)<=TIME_m & TIME_m<=EPOCHS(ii,2) ))) ...
                    - nanmean(PSTH_wrng_M(find( EPOCHS(ii,1)<=TIME_m & TIME_m<=EPOCHS(ii,2) ))) );
            end
        end
        
        
        
        if DEL_epoch<=2
            DIFF_actual = abs( nanmean(PSTH_corr_T(find( DELTA.START(DEL_epoch)<=TIME_t & TIME_t<=DELTA.END(DEL_epoch) ))) ...
                - nanmean(PSTH_wrng_T(find( DELTA.START(DEL_epoch)<=TIME_t & TIME_t<=DELTA.END(DEL_epoch) ))) );
        end
        if DEL_epoch>=3
            DIFF_actual = abs( nanmean(PSTH_corr_M(find( DELTA.START(DEL_epoch)<=TIME_m & TIME_m<=DELTA.END(DEL_epoch) ))) ...
                - nanmean(PSTH_wrng_M(find( DELTA.START(DEL_epoch)<=TIME_m & TIME_m<=DELTA.END(DEL_epoch) ))) );
        end
        
        
        
        
        subplot(4,4,4*(ite-1)+3)
        hist(DIFF); hold on;
        plot([DIFF_actual DIFF_actual],ylim,'-k');
        
        
        MEAN(ite,DEL_epoch) = nanmean(DIFF);
        MED(ite,DEL_epoch) = nanmedian(DIFF);
        STD(ite,DEL_epoch) = nanstd(DIFF);
        ACTUAL(ite,DEL_epoch) = DIFF_actual;
        
    end
    
    cd(Results_dir)
    filename = strcat('DELTA_ONLY_COGNITIVE_ERROR_',num2str(DEL_epoch));
    print(F, '-dpdf', filename, '-r400')
    
    
end



DELTA_stats.DIFF.MEAN = MEAN;
DELTA_stats.DIFF.MED = MED;
DELTA_stats.DIFF.STD = STD;
DELTA_stats.DIFF.ACTUAL = ACTUAL;


save(POP_file,'DELTA_stats','-append');
save(MERGE_file,'DELTA_stats','-append');
save(ALLCELLS_file,'DELTA_stats','-append');


