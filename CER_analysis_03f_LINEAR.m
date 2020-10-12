

%% function CER_analysis_03f_LINEAR

% Calculates the Hand Preference Index
% written by naveen at cumc on 8/6/17


% HANDS seperated -------------
IND_L = find(Infos(1:CHANGE,9)==1);
IND_R = find(Infos(1:CHANGE,9)==0);
SIGNAL = Spikes.S;    
Colour1 = [0.2 0.2 0.2];
Colour2 = [0.7 0.7 0.7];
sigma = 20;


Align_code = 4;
Start_time = -700;
End_time = 800;
P_Left_T  = PSTHe_n(SIGNAL(IND_L,1),Infos(IND_L,Align_code),Start_time,End_time,sigma,Colour1,1,1);
P_Right_T = PSTHe_n(SIGNAL(IND_R,1),Infos(IND_R,Align_code),Start_time,End_time,sigma,Colour2,1,1);


Align_code = 11;
Start_time = -500;
End_time = 700;
P_Left_M  = PSTHe_n(SIGNAL(IND_L,1),Infos(IND_L,Align_code),Start_time,End_time,sigma,Colour1,1,1);
P_Right_M = PSTHe_n(SIGNAL(IND_R,1),Infos(IND_R,Align_code),Start_time,End_time,sigma,Colour2,1,1);

% % % % Normalizing the neural response---------
% % % MIN = nanmin([nanmin(P_Left_T(1,:)), nanmin(P_Right_T(1,:)), nanmin(P_Left_M(1,:)), nanmin(P_Right_M(1,:)) ]);
% % % P_LEFT_T=P_Left_T(1,:)-MIN;
% % % P_RIGHT_T=P_Right_T(1,:)-MIN;
% % % P_LEFT_M=P_Left_M(1,:)-MIN;
% % % P_RIGHT_M=P_Right_M(1,:)-MIN;
% % % 
% % % MAX = nanmax([nanmax(P_Left_T(1,:)), nanmax(P_Right_T(1,:)), nanmax(P_Left_M(1,:)), nanmax(P_Right_M(1,:)) ]);
% % % 
% % % P_LEFT_T=P_LEFT_T/MAX;
% % % P_RIGHT_T=P_RIGHT_T/MAX;
% % % P_LEFT_M=P_LEFT_M/MAX;
% % % P_RIGHT_M=P_RIGHT_M/MAX;
% % % 
% % % % Then taking the difference---------
% % % DIFF_T = nanmean(P_LEFT_T-P_RIGHT_T)/nanmean(P_LEFT_T+P_RIGHT_T);
% % % DIFF_M = nanmean(P_LEFT_M-P_RIGHT_M)/nanmean(P_LEFT_M+P_RIGHT_M);
% % % 
% % % HPI = nanmean([DIFF_T,DIFF_M]);


DIFF_T = nanmean(P_Left_T(1,:)-P_Right_T(1,:))/nanmean(P_Left_T(1,:)+P_Right_T(1,:));
DIFF_M = nanmean(P_Left_M(1,:)-P_Right_M(1,:))/nanmean(P_Left_M(1,:)+P_Right_M(1,:));

HPI = nanmean([DIFF_T,DIFF_M]); % -1:left, +1: right


save(MERGE_file,'HPI','-append');
save(ALLCELLS_file,'HPI','-append');
save(POP_file,'HPI','-append');









% % % %% basic raster ----------------------------
% % % 
% % % F = figure();
% % % 
% % % subplot(3,2,1)
% % % Raster_n(Spikes.S(1:CHANGE-5,1),Infos(1:CHANGE-5,4),-500,1000,[0.4 0.4 0.4],0.25); % Aligned to target
% % % xlim([-400 200])
% % % 
% % % subplot(3,2,3)
% % % P = PSTHe_n(Spikes.S(1:CHANGE-5,1),Infos(1:CHANGE-5,4),-500,1000,20,[0.2 0.2 0.2],1,0);
% % % xlim([-400 200])
% % % time = -500:1000;
% % % ind = find(time>=-200 & time<=400);
% % % MIN = nanmin(P(1,ind)-P(2,ind));
% % % MAX = nanmax(P(1,ind)+P(2,ind));
% % % ylim([MIN MAX])
% % % 
% % % 
% % % 
% % % subplot(3,2,2)
% % % Raster_n(Spikes.S(1:CHANGE-5,1),Infos(1:CHANGE-5,11),-800,700,[0.4 0.4 0.4],0.25); % Aligned to target
% % % xlim([-200 400])
% % % 
% % % subplot(3,2,4)
% % % P = PSTHe_n(Spikes.S(1:CHANGE-5,1),Infos(1:CHANGE-5,11),-800,700,20,[0.2 0.2 0.2],1,0);
% % % xlim([-200 400])
% % % time = -800:700;
% % % ind = find(time>=-200 & time<=400);
% % % MIN = nanmin(P(1,ind)-P(2,ind));
% % % ylim([MIN MAX])
% % % 
% % % 
% % % suptitle(strcat('BASIC-',NOME));
% % % 
% % % cd(Results_dir)
% % % filename = strcat('BASIC_',NOME);
% % % print(F, '-dpdf', filename, '-r400')














%%%%%%%% GLOBAL MAX AND GLOBAL MIN %%%%%%%%%%%%%%%%%%



Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;


P_1 = PSTH_n(Spikes.S,Infos(:,4),Start_T,End_T,30,[0.2 0.2 0.2],1,1);
P_2 = PSTH_n(Spikes.S,Infos(:,11),Start_M,End_M,30,[0.2 0.2 0.2],1,1);


P1 = P_1([100: length(P_1)-100]);
P2 = P_2([100: length(P_2)-100]);


GLOBAL.MIN = nanmin([P1 P2]);
GLOBAL.MAX = nanmax([P1 P2]);


save(POP_file,'GLOBAL','-append');
save(MERGE_file,'GLOBAL','-append');
save(ALLCELLS_file,'GLOBAL','-append');







% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%% HAND AND RT? %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % clear Peaks Peak_pos
% % % % %     TIME = Start_M:End_M;
% % % % %     IND = find(TIME==-300):find(TIME==200);
% % % % % RT_ALL = Infos(:,14);
% % % % % for i=1:size(Infos,1)
% % % % %     P_ALL(i,:) = PSTH_ONE_n(Spikes.S{i,1},Infos(i,11),Start_M,End_M,60,[0.2 0.2 0.2],1,1);
% % % % % 
% % % % %     [Peaks(1,i) Peak_pos(1,i)] = nanmax(P_ALL(i,IND));
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % figure;
% % % % % plotyy(1:size(Infos,1),smooth(Peaks),1:size(Infos,1),smooth(RT_ALL));
% % % % % hold on;
% % % % % plot([CHANGE CHANGE],ylim,'-k')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%% RT and CW %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% CORR_IND = find(Infos(CHANGE:CHANGE+19,10)==1);
% WRNG_IND = find(Infos(CHANGE:CHANGE+19,10)==2);
% MAIN_IND = CHANGE-19:CHANGE;
% 
% CW_DELTA.RT_CORR = nanmean(Infos(CORR_IND,14));
% CW_DELTA.RT_WRNG = nanmean(Infos(WRNG_IND,14));
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE SHAPE and CW %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CW_DELTA.SHAPE_CORR = cell(4,1);
% CW_DELTA.SHAPE_WRNG = cell(4,1);
% CW_DELTA.SHAPE_MAIN = cell(4,1);
% CW_DELTA.SHAPE_CORRELATION_C_W = NaN(4,1);
% CW_DELTA.SHAPE_P_C_W = NaN(4,1);
% CW_DELTA.SHAPE_CORRELATION_C_M = NaN(4,1);
% CW_DELTA.SHAPE_P_C_M = NaN(4,1);
% CW_DELTA.SHAPE_CORRELATION_W_M = NaN(4,1);
% CW_DELTA.SHAPE_P_W_M = NaN(4,1);
% 
% 
% 
% 
% for index=1:2
%     if DELTA.EpochFlag(index)==1
%         Time = Start_T:End_T;
%         clear INDS
%         INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
%         
%         
%         % CORR ----------------------------------------------
%         temp_c = (Spikes.S(CORR_IND,1));
%         for jjj=1:length(temp_c)
%             TEMP_C{jjj,1} = temp_c{jjj,1}-Infos(CORR_IND(jjj),4);
%             SPK_IND{jjj,1}=find(TEMP_C{jjj,1}>=DELTA.START(index) & TEMP_C{jjj,1}<=DELTA.END(index));
%             SHAPES_C{jjj,:} = Shapes.S{CORR_IND(jjj),1}(SPK_IND{jjj,1},:);
%         end
%         
%         count_cc=0;
%         for kkk=1:size(SHAPES_C,1)
%             if ~isempty(SHAPES_C{kkk,1})
%                 for kkkk=1:size(SHAPES_C{kkk,1},1)
%                     count_cc=count_cc+1;
%                     SHAPES_C_fin(count_cc,:)=SHAPES_C{kkk,1}(kkkk,:);
%                 end
%             end
%         end
%         
%         % WRNG ----------------------------------------------
%         temp_w = (Spikes.S(WRNG_IND,1));
%         for jjj=1:length(temp_w)
%             TEMP_W{jjj,1} = temp_w{jjj,1}-Infos(WRNG_IND(jjj),4);
%             SPK_IND{jjj,1}=find(TEMP_W{jjj,1}>=DELTA.START(index) & TEMP_W{jjj,1}<=DELTA.END(index));
%             SHAPES_W{jjj,:} = Shapes.S{WRNG_IND(jjj),1}(SPK_IND{jjj,1},:);
%         end
%         
%         count_ww=0;
%         for kkk=1:size(SHAPES_W,1)
%             if ~isempty(SHAPES_W{kkk,1})
%                 for kkkk=1:size(SHAPES_W{kkk,1},1)
%                     count_ww=count_ww+1;
%                     SHAPES_W_fin(count_ww,:)=SHAPES_W{kkk,1}(kkkk,:);
%                 end
%             end
%         end
%         
%         
%          % MAIN ----------------------------------------------
%         temp_m = (Spikes.S(MAIN_IND,1));
%         for jjj=1:length(temp_m)
%             TEMP_M{jjj,1} = temp_m{jjj,1}-Infos(MAIN_IND(jjj),4);
%             SPK_IND{jjj,1}=find(TEMP_M{jjj,1}>=DELTA.START(index) & TEMP_M{jjj,1}<=DELTA.END(index));
%             SHAPES_M{jjj,:} = Shapes.S{MAIN_IND(jjj),1}(SPK_IND{jjj,1},:);
%         end
%         
%         count_mm=0;
%         for kkk=1:size(SHAPES_M,1)
%             if ~isempty(SHAPES_M{kkk,1})
%                 for kkkk=1:size(SHAPES_M{kkk,1},1)
%                     count_mm=count_mm+1;
%                     SHAPES_M_fin(count_mm,:)=SHAPES_M{kkk,1}(kkkk,:);
%                 end
%             end
%         end
%         
%         CW_DELTA.SHAPE_CORR{index,:} = nanmean(SHAPES_C_fin);
%         CW_DELTA.SHAPE_WRNG{index,:} = nanmean(SHAPES_W_fin);
%         [CW_DELTA.SHAPE_CORRELATION_C_W(index,1), CW_DELTA.SHAPE_P_C_W(index,1)] = corr(nanmean(SHAPES_C_fin)',nanmean(SHAPES_W_fin)');
%         
%         CW_DELTA.SHAPE_MAIN{index,:} = nanmean(SHAPES_M_fin);
%         [CW_DELTA.SHAPE_CORRELATION_C_M(index,1), CW_DELTA.SHAPE_P_C_M(index,1)] = corr(nanmean(SHAPES_C_fin)',nanmean(SHAPES_M_fin)');
%         [CW_DELTA.SHAPE_CORRELATION_W_M(index,1), CW_DELTA.SHAPE_P_W_M(index,1)] = corr(nanmean(SHAPES_M_fin)',nanmean(SHAPES_W_fin)');
%         
%         
%     end
% end
% 
% 
% 
% 
% index=3;
% if DELTA.EpochFlag(index)==1
%     Time = Start_M:End_M;
%     clear INDS
%     INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
%     
%     % CORR ----------------------------------------------
%     temp_c = (Spikes.S(CORR_IND,1));
%     for jjj=1:length(temp_c)
%         TEMP_C{jjj,1} = temp_c{jjj,1}-Infos(CORR_IND(jjj),11);
%         SPK_IND{jjj,1}=find(TEMP_C{jjj,1}>=DELTA.START(index) & TEMP_C{jjj,1}<=DELTA.END(index));
%         SHAPES_C{jjj,:} = Shapes.S{CORR_IND(jjj),1}(SPK_IND{jjj,1},:);
%     end
%     
%     count_cc=0;
%     for kkk=1:size(SHAPES_C,1)
%         if ~isempty(SHAPES_C{kkk,1})
%             for kkkk=1:size(SHAPES_C{kkk,1},1)
%                 count_cc=count_cc+1;
%                 SHAPES_C_fin(count_cc,:)=SHAPES_C{kkk,1}(kkkk,:);
%             end
%         end
%     end
%     
%     % WRNG ----------------------------------------------
%     temp_w = (Spikes.S(WRNG_IND,1));
%     for jjj=1:length(temp_w)
%         TEMP_W{jjj,1} = temp_w{jjj,1}-Infos(WRNG_IND(jjj),11);
%         SPK_IND{jjj,1}=find(TEMP_W{jjj,1}>=DELTA.START(index) & TEMP_W{jjj,1}<=DELTA.END(index));
%         SHAPES_W{jjj,:} = Shapes.S{WRNG_IND(jjj),1}(SPK_IND{jjj,1},:);
%     end
%     
%     count_ww=0;
%     for kkk=1:size(SHAPES_W,1)
%         if ~isempty(SHAPES_W{kkk,1})
%             for kkkk=1:size(SHAPES_W{kkk,1},1)
%                 count_ww=count_ww+1;
%                 SHAPES_W_fin(count_ww,:)=SHAPES_W{kkk,1}(kkkk,:);
%             end
%         end
%     end
%     
%     
%     
%    % MAIN ----------------------------------------------
%         temp_m = (Spikes.S(MAIN_IND,1));
%         for jjj=1:length(temp_m)
%             TEMP_M{jjj,1} = temp_m{jjj,1}-Infos(MAIN_IND(jjj),11);
%             SPK_IND{jjj,1}=find(TEMP_M{jjj,1}>=DELTA.START(index) & TEMP_M{jjj,1}<=DELTA.END(index));
%             SHAPES_M{jjj,:} = Shapes.S{MAIN_IND(jjj),1}(SPK_IND{jjj,1},:);
%         end
%         
%         count_mm=0;
%         for kkk=1:size(SHAPES_M,1)
%             if ~isempty(SHAPES_M{kkk,1})
%                 for kkkk=1:size(SHAPES_M{kkk,1},1)
%                     count_mm=count_mm+1;
%                     SHAPES_M_fin(count_mm,:)=SHAPES_M{kkk,1}(kkkk,:);
%                 end
%             end
%         end
%         
%         CW_DELTA.SHAPE_CORR{index,:} = nanmean(SHAPES_C_fin);
%         CW_DELTA.SHAPE_WRNG{index,:} = nanmean(SHAPES_W_fin);
%         [CW_DELTA.SHAPE_CORRELATION_C_W(index,1), CW_DELTA.SHAPE_P_C_W(index,1)] = corr(nanmean(SHAPES_C_fin)',nanmean(SHAPES_W_fin)');
%         
%         CW_DELTA.SHAPE_MAIN{index,:} = nanmean(SHAPES_M_fin);
%         [CW_DELTA.SHAPE_CORRELATION_C_M(index,1), CW_DELTA.SHAPE_P_C_M(index,1)] = corr(nanmean(SHAPES_C_fin)',nanmean(SHAPES_M_fin)');
%         [CW_DELTA.SHAPE_CORRELATION_W_M(index,1), CW_DELTA.SHAPE_P_W_M(index,1)] = corr(nanmean(SHAPES_M_fin)',nanmean(SHAPES_W_fin)');
%         
% end
% 
% 
% index=4;
% if DELTA.EpochFlag(index)==1
%     Time = Start_M:End_M;
%     clear INDS
%     INDS = find(DELTA.START(index)<=Time & Time<=DELTA.END(index));
%     % CORR ----------------------------------------------
%     temp_c = (Spikes.S(CORR_IND,1));
%     for jjj=1:length(temp_c)
%         TEMP_C{jjj,1} = temp_c{jjj,1}-Infos(CORR_IND(jjj),11);
%         SPK_IND{jjj,1}=find(TEMP_C{jjj,1}>=DELTA.START(index) & TEMP_C{jjj,1}<=DELTA.END(index));
%         SHAPES_C{jjj,:} = Shapes.S{CORR_IND(jjj),1}(SPK_IND{jjj,1},:);
%     end
%     
%     count_cc=0;
%     for kkk=1:size(SHAPES_C,1)
%         if ~isempty(SHAPES_C{kkk,1})
%             for kkkk=1:size(SHAPES_C{kkk,1},1)
%                 count_cc=count_cc+1;
%                 SHAPES_C_fin(count_cc,:)=SHAPES_C{kkk,1}(kkkk,:);
%             end
%         end
%     end
%     
%     % WRNG ----------------------------------------------
%     temp_w = (Spikes.S(WRNG_IND,1));
%     for jjj=1:length(temp_w)
%         TEMP_W{jjj,1} = temp_w{jjj,1}-Infos(WRNG_IND(jjj),11);
%         SPK_IND{jjj,1}=find(TEMP_W{jjj,1}>=DELTA.START(index) & TEMP_W{jjj,1}<=DELTA.END(index));
%         SHAPES_W{jjj,:} = Shapes.S{WRNG_IND(jjj),1}(SPK_IND{jjj,1},:);
%     end
%     
%     count_ww=0;
%     for kkk=1:size(SHAPES_W,1)
%         if ~isempty(SHAPES_W{kkk,1})
%             for kkkk=1:size(SHAPES_W{kkk,1},1)
%                 count_ww=count_ww+1;
%                 SHAPES_W_fin(count_ww,:)=SHAPES_W{kkk,1}(kkkk,:);
%             end
%         end
%     end
%     
%    
%     
%    % MAIN ----------------------------------------------
%         temp_m = (Spikes.S(MAIN_IND,1));
%         for jjj=1:length(temp_m)
%             TEMP_M{jjj,1} = temp_m{jjj,1}-Infos(MAIN_IND(jjj),11);
%             SPK_IND{jjj,1}=find(TEMP_M{jjj,1}>=DELTA.START(index) & TEMP_M{jjj,1}<=DELTA.END(index));
%             SHAPES_M{jjj,:} = Shapes.S{MAIN_IND(jjj),1}(SPK_IND{jjj,1},:);
%         end
%         
%         count_mm=0;
%         for kkk=1:size(SHAPES_M,1)
%             if ~isempty(SHAPES_M{kkk,1})
%                 for kkkk=1:size(SHAPES_M{kkk,1},1)
%                     count_mm=count_mm+1;
%                     SHAPES_M_fin(count_mm,:)=SHAPES_M{kkk,1}(kkkk,:);
%                 end
%             end
%         end
%         
%         CW_DELTA.SHAPE_CORR{index,:} = nanmean(SHAPES_C_fin);
%         CW_DELTA.SHAPE_WRNG{index,:} = nanmean(SHAPES_W_fin);
%         [CW_DELTA.SHAPE_CORRELATION_C_W(index,1), CW_DELTA.SHAPE_P_C_W(index,1)] = corr(nanmean(SHAPES_C_fin)',nanmean(SHAPES_W_fin)');
%         
%         CW_DELTA.SHAPE_MAIN{index,:} = nanmean(SHAPES_M_fin);
%         [CW_DELTA.SHAPE_CORRELATION_C_M(index,1), CW_DELTA.SHAPE_P_C_M(index,1)] = corr(nanmean(SHAPES_C_fin)',nanmean(SHAPES_M_fin)');
%         [CW_DELTA.SHAPE_CORRELATION_W_M(index,1), CW_DELTA.SHAPE_P_W_M(index,1)] = corr(nanmean(SHAPES_M_fin)',nanmean(SHAPES_W_fin)');
%         
% end
% 
% 
% 
% save(POP_file,'CW_DELTA','-append');
% save(MERGE_file,'CW_DELTA','-append');
% save(ALLCELLS_file,'CW_DELTA','-append');
% 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIMULUS PREFERENCE ------------------------------------
% 
% % if ~strcmp(NOME(1),'P')
%     % STEP 1: FIND THE TWO DIFF STIMULI
%     
% %     SYM = unique(Infos(CHANGE:CHANGE+10,5));
% SYM = [999; 888];
%     
%     % SYMBOLS seperated -------------
%     IND_ONE = find(Infos(1:CHANGE,5)==SYM(1));
%     IND_TWO = find(Infos(1:CHANGE,5)==SYM(2));
%     
% % end
% % 
% % if strcmp(NOME(1),'P')
% %   
% % temp1 = NaN(size(Infos,1),1);
% % Stim  = NaN(size(Infos,1),1);
% % flag  = zeros(4,1);
% % temp2 = 0;
% % temp3 = NaN(size(Infos,1),1);
% % temps = NaN(4,1);
% % 
% % for i=1:size(Infos,1)
% % 
% %     if Infos(i,5)==160
% %         temp1(i) = (4*(Infos(i,6))) + (9*(Infos(i,7))) + (16*(Infos(i,8)));
% %     end
% % 
% %     if Infos(i,5)==160
% %         if (temp1(1) == temp1(i))
% %             Stim(i)=1; else Stim(i)=2; end
% %     end
% % 
% % 
% %     if (Infos(i,5)>160) && (flag(1)==0)
% %         flag(1)=1;
% %         temps(1) = Infos(i,5);
% %     end
% % 
% %     if (Infos(i,5)>160) && (flag(1)==1)
% %         temp3(i) = Infos(i,5);
% % 
% %         if (temp3(i) ~= temps(1)) && (flag(2)==0)
% %             temps(2)=temp3(i);
% %             flag(2)=1;
% %         end
% % 
% %         if (temp3(i) ~= temps(1)) && (temp3(i) ~= temps(2)) && (flag(2)==1) && (flag(3)==0)
% %             temps(3)=temp3(i);
% %             flag(3)=1;
% %         end
% % 
% %         if (temp3(i) ~= temps(1)) && (temp3(i) ~= temps(2)) && (temp3(i) ~= temps(3)) && (flag(3)==1) && (flag(4)==0)
% %             temps(4)=temp3(i);
% %             flag(4)=1;
% %         end
% % 
% %     end
% % 
% % 
% %     if (Infos(i,5)>160) && (flag(1)==1)
% %         if (temp3(i) == temps(1))
% %             Stim(i)=3;
% %         elseif (temp3(i) == temps(2))
% %             Stim(i)=4;
% %         elseif (temp3(i) == temps(3))
% %             Stim(i)=5;
% %         elseif (temp3(i) == temps(4))
% %             Stim(i)=6;
% %         end
% %     end
% % 
% % 
% % end  
% %     
% % end
% %  
% 
% 
% 
% % 
% % SIGNAL = Spikes.S;    
% % Colour1 = [0.2 0.2 0.2];
% % Colour2 = [0.7 0.7 0.7];
% % sigma = 20;
% % 
% % 
% % Align_code = 4;
% % Start_time = -700;
% % End_time = 800;
% % P_ONE_T  = PSTH_n(SIGNAL(IND_ONE,1),Infos(IND_ONE,Align_code),Start_time,End_time,sigma,Colour1,1,1);
% % P_TWO_T = PSTH_n(SIGNAL(IND_TWO,1),Infos(IND_TWO,Align_code),Start_time,End_time,sigma,Colour2,1,1);
% % 
% % 
% % Align_code = 11;
% % Start_time = -500;
% % End_time = 700;
% % P_ONE_M  = PSTH_n(SIGNAL(IND_ONE,1),Infos(IND_ONE,Align_code),Start_time,End_time,sigma,Colour1,1,1);
% % P_TWO_M = PSTH_n(SIGNAL(IND_TWO,1),Infos(IND_TWO,Align_code),Start_time,End_time,sigma,Colour2,1,1);
% % 
% % 
% % 
% % diff_T = nanmean(P_ONE_T(1,:)-P_TWO_T(1,:))/nanmean(P_ONE_T(1,:)+P_TWO_T(1,:));
% % diff_M = nanmean(P_ONE_M(1,:)-P_TWO_M(1,:))/nanmean(P_ONE_M(1,:)+P_TWO_M(1,:));
% % 
% % TPI = nanmean([diff_T,diff_M]); % -1:left, +1: right
% % 
% % 
% % save(MERGE_file,'TPI','-append');
% % save(ALLCELLS_file,'TPI','-append');
% % save(POP_file,'TPI','-append');
% % 
% % % end
% % 
% 
% 
% 
% 
% 
% 
% 





% % % %% %%%%%%%%%%% TRIAL BASED LEARNING
% % % 
% % % outcome = Infos(CHANGE:LEARNT,10);
% % % LEN = length(outcome);
% % % LEN = round(0.5*LEN);
% % % outcome = outcome(1:LEN);
% % % outcome(find(outcome==2))=0;
% % % aftW = find(outcome(1:end-1)==0)+1;
% % % aftC = find(outcome(1:end-1)==1)+1;
% % % 
% % % AFT.W = nanmean(outcome(aftW));
% % % AFT.C = nanmean(outcome(aftC));
% % % 
% % % 
% % % save(MERGE_file,'AFT','-append');
% % % save(ALLCELLS_file,'AFT','-append');
% % % save(POP_file,'AFT','-append');




