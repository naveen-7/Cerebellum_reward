

% RT_SUPP
% code to analyse the firing rate changes between correct trials when there
% was no change in RT

% written by naveen at cumc on 9/25/17


disp('!!!!!  CER_analysis_POPULATION_LINEAR_n_RT_SUPP has started running  !!!!!')


% Setup directories--------------------------------------------------------

codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data\POP_CELLS';

cd(data_dir)
mkdir('Results');
Results_pop = strcat(data_dir,'/RESULTS');
PATHNAME    = strcat(data_dir,'/');


clear FILES_RT;

RT_CELLS = {'POP_sl080817_01N';
'POP_sl072617_01N';
'POP_sl071417_03N';
'POP_sl071417_01N';
'POP_sl071117_04N';
'POP_sl071117_01N';
'POP_sl070617_01N';
'POP_sl063017_01N';
'POP_pr080315_001';
'POP_pr070115_005';
'POP_pr062214_002';
'POP_pr062214_001';
'POP_br121916_02N';
'POP_br102716_02R';
'POP_br052416_03N';
'POP_br040416_005';
'POP_br031916_002'};



count = 0;
for kk=1:length(RT_CELLS)
    str = cell2mat(RT_CELLS(kk,:));
    if (length(str)>5)
        if (strcmp(str(1:4),'POP_'))
            count = count+1;
            FILES_RT(count,1:length(RT_CELLS(kk,:) )) = RT_CELLS(kk,:);
        end
    end
end



for kk=1:length(FILES_RT)
    str = cell2mat(FILES_RT(kk,:));
    temp = str(15:length(str)-4);
end



RT_SUPP_POP      = cell(size(FILES_RT,1),1);
RT_GLOBAL_POP    = cell(size(FILES_RT,1),1);
RT_actual_POP    = cell(size(FILES_RT,1),1);
LC_actual_POP    = cell(size(FILES_RT,1),1);
CHANGE_POP       = NaN(size(FILES_RT,1),1);

ite=0;
for i=1 : size(FILES_RT,1)
    i
    clear HPI RT TPI GLOBAL NOME TASK_TYPE REWARD CW_DELTA DELTA STATE VARIANCE LEARNING;
    
    Curr_file = cell2mat(strcat(PATHNAME,FILES_RT(i,:)))
    Curr_file_Name = cell2mat(FILES_RT(i,:));
    load(Curr_file);
    
    
    tempname  = FILES_RT(i,:);
    tempname = tempname{1,1};
    tempname = strcat('DATA_',(tempname(5:length(tempname))));
    load(fullfile('C:','NAVEEN_Work','Cerebellum','Data','ALL_CELLS',tempname),'Infos','CHANGE');
    
    ite=ite+1;
    
    CELL_NOME{ite,1}     = NOME;  %subtype
    TASK_TYPE_POP{ite,1} = TASK_TYPE;
    
    if exist('LEARNING')
        LEARNING_POP{ite,1} = LEARNING;
        if strcmp(LEARNING,'N') DELTA=[]; STATE=[]; VARIANCE=[];   end
    end
    
    if exist('CW_DELTA')   RT_SUPP_POP{ite,1} = CW_DELTA;     end
    if exist('GLOBAL')     RT_GLOBAL_POP{ite,1} = GLOBAL;     end
    RT_actual_POP{ite,1} = Infos(:,14);
    LC_actual_POP{ite,1} = Infos(:,10);
    CHANGE_POP(ite,1) = CHANGE;
    
    clear Infos CHANGE;
 
end




disp('*************************************************')
disp('*************************************************')
disp('ALL FILES WERE LOADED SUCCESSFULLY')






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% WITH MAIN STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INITIALIZATION ------------------------------------------------
MT_POP_NORM_RT = cell(size(FILES_RT,1),1);

% EXTRACTION ----------------------------------------------------
for i=1 : size(FILES_RT,1)
    if  ~isempty(RT_SUPP_POP{i, 1} )
        if  isfield(RT_SUPP_POP{i,1},'statsMT')
            MT_POP_NORM_RT{i,1} = RT_SUPP_POP{i, 1}.statsMT;
            MT_POP_NORM_RTstd{i,1} = RT_SUPP_POP{i, 1}.statsMTstd;
        end
    end
end

% 
% 
% for i=1:size(FILES_RT,1)
%     if  ~isempty(MT_POP_NORM_RT{i,1})
%         MT_POP_NORM_RT{i,1} = MT_POP_NORM_RT{i,1}-RT_GLOBAL_POP{i,1}.MIN;
%         MT_POP_NORM_RT{i,1}(MT_POP_NORM_RT{i,1}<0)=0; 
%         MT_POP_NORM_RT{i,1} = MT_POP_NORM_RT{i,1}/RT_GLOBAL_POP{i,1}.MAX;
%     end
% end
% MIN=0; MAX=1;


% % % 
% % % for i=1:size(FILES_RT,1)
% % %     if  ~isempty(MT_POP_NORM{i,1})
% % %         MT_POP_NORM{i,1} = MT_POP_NORM{i,1}/nanmean([GLOBAL_POP{i,1}.MAX GLOBAL_POP{i,1}.MIN]);
% % %     end
% % % end
% % % MIN=0; MAX=3;






% INITIALIZATION ------------------------------------------------
CN_POP_NORM_RT = cell(size(FILES_RT,1),1);
WN_POP_NORM_RT = cell(size(FILES_RT,1),1);
NC_POP_NORM_RT = cell(size(FILES_RT,1),1);
NW_POP_NORM_RT = cell(size(FILES_RT,1),1);

WRNG_NORM_RT = NaN(size(FILES_RT,1),4);
CORR_NORM_RT = NaN(size(FILES_RT,1),4);
P_NORM_RT = NaN(size(FILES_RT,1),4);

P_POP_NORM_RT = cell(size(FILES_RT,1),1);


% EXTRACTION ----------------------------------------------------
for i=1 : size(FILES_RT,1)
    if  ~isempty(RT_SUPP_POP{i, 1} )
        if  isfield(RT_SUPP_POP{i,1},'statsCN')
            CN_POP_NORM_RT{i,1} = RT_SUPP_POP{i, 1}.statsCN;
            WN_POP_NORM_RT{i,1} = RT_SUPP_POP{i, 1}.statsWN;
            NC_POP_NORM_RT{i,1} = RT_SUPP_POP{i, 1}.statsNC;
            NW_POP_NORM_RT{i,1} = RT_SUPP_POP{i, 1}.statsNW;
            
            CN_POP_NORM_RTstd{i,1} = RT_SUPP_POP{i, 1}.statsCNstd;
            WN_POP_NORM_RTstd{i,1} = RT_SUPP_POP{i, 1}.statsWNstd;
            NC_POP_NORM_RTstd{i,1} = RT_SUPP_POP{i, 1}.statsNCstd;
            NW_POP_NORM_RTstd{i,1} = RT_SUPP_POP{i, 1}.statsNWstd;
            
            P_POP_NORM_RT{i,1} = RT_SUPP_POP{i, 1}.statsP;
        end
    end
end

% 
% %%% normalizing min to max
% for i=1:size(FILES_RT,1)
%     if  ~isempty(WN_POP_NORM_RT{i,1})
%         WN_POP_NORM_RT{i,1} = WN_POP_NORM_RT{i,1}-RT_GLOBAL_POP{i,1}.MIN;
%         WN_POP_NORM_RT{i,1}(WN_POP_NORM_RT{i,1}<0)=0; 
%         WN_POP_NORM_RT{i,1} = WN_POP_NORM_RT{i,1}/RT_GLOBAL_POP{i,1}.MAX;
%         
%         NW_POP_NORM_RT{i,1} = NW_POP_NORM_RT{i,1}-RT_GLOBAL_POP{i,1}.MIN;
%         NW_POP_NORM_RT{i,1}(NW_POP_NORM_RT{i,1}<0)=0; 
%         NW_POP_NORM_RT{i,1} = NW_POP_NORM_RT{i,1}/RT_GLOBAL_POP{i,1}.MAX;
%         
%         CN_POP_NORM_RT{i,1} = CN_POP_NORM_RT{i,1}-RT_GLOBAL_POP{i,1}.MIN;
%         CN_POP_NORM_RT{i,1}(CN_POP_NORM_RT{i,1}<0)=0; 
%         CN_POP_NORM_RT{i,1} = CN_POP_NORM_RT{i,1}/RT_GLOBAL_POP{i,1}.MAX;
%         
%         NC_POP_NORM_RT{i,1} = NC_POP_NORM_RT{i,1}-RT_GLOBAL_POP{i,1}.MIN;
%         NC_POP_NORM_RT{i,1}(NC_POP_NORM_RT{i,1}<0)=0; 
%         NC_POP_NORM_RT{i,1} = NC_POP_NORM_RT{i,1}/RT_GLOBAL_POP{i,1}.MAX;
%     end
% end
% MIN = 0; MAX = 1;
% 












%%% JUST FOR THOSE WITH MAIN

VAL=0;
count_WITHMAIN_RT=0;
clear WRNG_PREF_WITHMAIN_RT CORR_PREF_WITHMAIN_RT WITHMAIN_RT WITHMAIN_RTstd
clear WRNG_PREF_WITHMAIN_C_RT WRNG_PREF_WITHMAIN_W_RT CORR_PREF_WITHMAIN_C_RT CORR_PREF_WITHMAIN_W_RT
for i=1:size(FILES_RT,1)
    if  ~isempty(MT_POP_NORM_RT{i,1})
        i
        WRNG_WITHMAIN_RT(i,1) = WN_POP_NORM_RT{i,1}(1);
        WRNG_WITHMAIN_RT(i,2) = WN_POP_NORM_RT{i,1}(2);
        WRNG_WITHMAIN_RT(i,3) = WN_POP_NORM_RT{i,1}(3);
        WRNG_WITHMAIN_RT(i,4) = NW_POP_NORM_RT{i,1}(4);
        
        CORR_WITHMAIN_RT(i,1) = CN_POP_NORM_RT{i,1}(1);
        CORR_WITHMAIN_RT(i,2) = CN_POP_NORM_RT{i,1}(2);
        CORR_WITHMAIN_RT(i,3) = CN_POP_NORM_RT{i,1}(3);
        CORR_WITHMAIN_RT(i,4) = NC_POP_NORM_RT{i,1}(4);
        
        MAIN_RT(i,1) = MT_POP_NORM_RT{i,1}(1);
        MAIN_RT(i,2) = MT_POP_NORM_RT{i,1}(2);
        MAIN_RT(i,3) = MT_POP_NORM_RT{i,1}(3);
        MAIN_RT(i,4) = MT_POP_NORM_RT{i,1}(4);
        
        
        CORR_WITHMAIN_RTstd(i,1) = CN_POP_NORM_RTstd{i,1}(1);
        CORR_WITHMAIN_RTstd(i,2) = CN_POP_NORM_RTstd{i,1}(2);
        CORR_WITHMAIN_RTstd(i,3) = CN_POP_NORM_RTstd{i,1}(3);
        CORR_WITHMAIN_RTstd(i,4) = NC_POP_NORM_RTstd{i,1}(4);
        
        MAIN_RTstd(i,1) = MT_POP_NORM_RTstd{i,1}(1);
        MAIN_RTstd(i,2) = MT_POP_NORM_RTstd{i,1}(2);
        MAIN_RTstd(i,3) = MT_POP_NORM_RTstd{i,1}(3);
        MAIN_RTstd(i,4) = MT_POP_NORM_RTstd{i,1}(4);
       
        %%% JUST FOR THOSE WITH MAIN
        
         for j=1:4
            if ~isnan(WRNG_WITHMAIN_RT(i,j)) & ~isnan(CORR_WITHMAIN_RT(i,j)) & ~isnan(MAIN_RT(i,j))
                 if MAIN_RT(i,j)>CORR_WITHMAIN_RT(i,j)+VAL & MAIN_RT(i,j)<WRNG_WITHMAIN_RT(i,j)-VAL | ...
                         MAIN_RT(i,j)<CORR_WITHMAIN_RT(i,j)+VAL & MAIN_RT(i,j)>WRNG_WITHMAIN_RT(i,j)-VAL
                    
                   count_WITHMAIN_RT=count_WITHMAIN_RT+1;
                   WITHMAIN_RT{count_WITHMAIN_RT,1} = [MAIN_RT(i,j) CORR_WITHMAIN_RT(i,j)];
                   WITHMAIN_RTstd{count_WITHMAIN_RT,1} = [MAIN_RTstd(i,j) CORR_WITHMAIN_RTstd(i,j)];
                    
                 end
            end
        end    
    end
end

count_WITHMAIN_RT


clear TTEST
WITHMAIN_RT_mat = cell2mat(WITHMAIN_RT);
WITHMAIN_RTstd_mat = cell2mat(WITHMAIN_RTstd);

for iii = 1:size(WITHMAIN_RT_mat,1)
    TTEST_NOCHANGE(iii,1) = ttest_n(WITHMAIN_RT_mat(iii,1),WITHMAIN_RT_mat(iii,2),WITHMAIN_RTstd_mat(iii,1),WITHMAIN_RTstd_mat(iii,2),20,20);
end

% % % % TTEST_NOCHANGE = TTEST_NOCHANGE(~TTEST_NOCHANGE==0);
% % % 
% % % TTEST_NOCHANGE(TTEST_NOCHANGE==0)=10.^-randi([20 50],length(TTEST_NOCHANGE(TTEST_NOCHANGE==0)),1);
% % % TTEST_NOCHANGE(TTEST_NOCHANGE>0.01)=NaN;
% % % TTEST_NOCHANGE = TTEST_NOCHANGE(~isnan(TTEST_NOCHANGE));

%%%% DOING RT for nochange -------

for i=1:size(RT_actual_POP,1)
    temp_rt = RT_actual_POP{i,1};
    RT_toplot(i,1)= nanmean(temp_rt(CHANGE_POP(i,1)-10:CHANGE_POP(i,1)));
    RT_toplot(i,2)= nanmean(temp_rt(CHANGE_POP(i,1)+1:CHANGE_POP(i,1)+11));  
    clear temp_rt;
end


NOchange_RT(1,1) = nanmean(RT_toplot(:,1));
NOchange_RT(1,2) = nanmean(RT_toplot(:,2));

NOchange_RT(2,1) = nanstd(RT_toplot(:,1))/sqrt(size(RT_toplot,1));
NOchange_RT(2,2) = nanstd(RT_toplot(:,2))/sqrt(size(RT_toplot,1));

Nochange_stats_RT = ttest_NN(RT_toplot(:,1),RT_toplot(:,2));




%%%% DOING LC for nochange -------

for i=1:size(LC_actual_POP,1)
    temp_lc = LC_actual_POP{i,1};
    temp_lc(temp_lc==2)=0;
    LC_toplot(i,1)= nanmean(temp_lc(CHANGE_POP(i,1)-10:CHANGE_POP(i,1)))*100;
    LC_toplot(i,2)= nanmean(temp_lc(CHANGE_POP(i,1)+1:CHANGE_POP(i,1)+11))*100;  
    clear temp_lc;
end


NOchange_LC(1,1) = nanmean(LC_toplot(:,1));
NOchange_LC(1,2) = nanmean(LC_toplot(:,2));

NOchange_LC(2,1) = nanstd(LC_toplot(:,1))/sqrt(size(LC_toplot,1));
NOchange_LC(2,2) = nanstd(LC_toplot(:,2))/sqrt(size(LC_toplot,1));


Nochange_stats_LC = ttest_NN(LC_toplot(:,1),LC_toplot(:,2));







RT_COLOUR = [255 102 102]/255;




% F = figure
subplot(4,4,4)
hold on;
b = bar(1,NOchange_LC(1,1))
b.FaceColor = [0.5 0.5 0.5];
b = bar(2,NOchange_LC(1,2))
b.FaceColor = [0.5 0.5 0.5];

errorbar([1 2], ...
    [NOchange_LC(1,1) NOchange_LC(1,2)], ...
    [NOchange_LC(2,1) NOchange_LC(2,2)], ...
    '.K')
box off;
xlim([0 3])

MIN = 50; MAX = 100;
plot([1.2 1.8],[MAX-5 MAX-5],'-K');
text(1.2, MAX,star_n(Nochange_stats_LC));
ylim([MIN MAX]);
ylabel('% correct (ms)');
set(gca,'xTick',[1 2])
set(gca, 'XTickLabel', {'Bef','Aft'},'fontsize',7);
set(gca,'linewidth',0.6,'fontsize',8,'fontweight','normal')



subplot(4,4,8)
hold on;
b = bar(1,NOchange_RT(1,1))
b.FaceColor = RT_COLOUR;
b = bar(2,NOchange_RT(1,2))
b.FaceColor = RT_COLOUR;

errorbar([1 2], ...
    [NOchange_RT(1,1) NOchange_RT(1,2)], ...
    [NOchange_RT(2,1) NOchange_RT(2,2)], ...
    '.K')
box off;
xlim([0 3])

MAX = 850; MIN=500;
plot([1.2 1.8],[MAX MAX],'-K');
text(1.2, MAX+50,star_n(Nochange_stats_RT));
ylim([MIN 1000]);
ylabel('RT (ms)');
set(gca,'xTick',[1 2])
set(gca, 'XTickLabel', {'Bef','Aft'},'fontsize',7);
set(gca,'linewidth',0.6,'fontsize',8,'fontweight','normal')



subplot(4,4,12)
hold on;

h = histogram(log10(TTEST_NOCHANGE),10);
% h = histogram((TTEST_YESCHANGE),20);
hold on;

h.FaceColor = RT_COLOUR;
box off;
xlim([-70 0])
ylim([0 8])
plot([-3 -3],ylim,'--k')
ylabel('log(p)')


% beeswarm_plot_n({log10(TTEST_NOCHANGE)},'color',RT_COLOUR,'stats','none','MS',20);
% plot(1,nanmean(log10(TTEST_NOCHANGE)),'ok');


% 
% filename = 'RT_SUPP_EXTRA';
% cd(Results_pop);
% print(F, '-dpdf', filename, '-r400');
% 
% clear F;









































