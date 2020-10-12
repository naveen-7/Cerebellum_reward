

% function CER_analysis_03e_LINEAR

% % % 
% % % 
% % % for index=1:4
% % %     if DELTA.EpochFlag(index)==1
% % %         
% % %         index
% % %         
% % %         if index==1 EPOCH_TYPE = 'B'; end
% % %         if index==2 EPOCH_TYPE = 'T'; end
% % %         if index==3 EPOCH_TYPE = 'M'; end
% % %         if index==4 EPOCH_TYPE = 'P'; end
% % %         
% % %         
% % %         EPOCH_START = DELTA.START(index);
% % %         EPOCH_END = DELTA.END(index);
% % %         
% % %         
% % %         
% % %         
% % %         % FIRST THINGS ------------
% % %         
% % %         % Start = EPOCH_START;
% % %         % End = EPOCH_END;
% % %         
% % %         if strcmp(EPOCH_TYPE,'P') | strcmp(EPOCH_TYPE,'M')
% % %             Align_code = 11;
% % %             ALIGN = 'M';
% % %         end
% % %         
% % %         if strcmp(EPOCH_TYPE,'B') | strcmp(EPOCH_TYPE,'T')
% % %             Align_code = 4;
% % %             ALIGN = 'T';
% % %         end
% % %         
% % %         if Align_code == 11
% % %             Start = -750; End = 750;
% % %         end
% % %         
% % %         if Align_code == 4
% % %             Start = -450; End = 1050;
% % %         end
% % %         
% % %         TIME = Start:End;
% % %         
% % %         
% % %         clear PSTH1 Matrix1 PSTH2 Matrix2 MAT
% % %         for i = 1:size(Infos,1)  % Running for all trials in each case
% % %             Sig = Spikes.S{i,1};
% % %             PSTH1 = PSTH_ONE_n(Sig,Infos(i,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
% % %             Matrix1(i,:) = PSTH1(find(TIME==EPOCH_START):find(TIME==EPOCH_END));
% % %             MAT(i,:) = PSTH1;
% % %             %     PSTH2 = PSTH_ONE_n(Sig,Infos(i,Align_code),Start_RANDOM,End_RANDOM,60,[0.9804    0.5020    0.4471]);
% % %             %     Matrix2(i,:) = PSTH2(find(TIME==Start_RANDOM):find(TIME==End_RANDOM));
% % %         end
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         %% Variance PRELIM -----------
% % %         %
% % %         % clear PSTH1 MAT
% % %         % for i = 1:size(Infos,1)  % Running for all trials in each case
% % %         %     Sig = Spikes.S{i,1};
% % %         %     PSTH1 = PSTH_ONE_n(Sig,Infos(i,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
% % %         %     MAT(i,:) = PSTH1;
% % %         % end
% % %         
% % %         BEF(1,:) = nanmean(MAT(1:CHANGE,:));
% % %         BEF(2,:) = nanstd(MAT(1:CHANGE,:));
% % %         BEF(3,:) = nanvar(MAT(1:CHANGE,:))./nanmean(MAT(1:CHANGE,:));
% % %         
% % %         DUR(1,:) = nanmean(MAT(CHANGE:LEARNT,:));
% % %         DUR(2,:) = nanstd(MAT(CHANGE:LEARNT,:));
% % %         DUR(3,:) = nanvar(MAT(CHANGE:LEARNT,:))./nanmean(MAT(CHANGE:LEARNT,:));
% % %         
% % %         AFT(1,:) = nanmean(MAT(LEARNT:size(Infos,1),:));
% % %         AFT(2,:) = nanstd(MAT(LEARNT:size(Infos,1),:));
% % %         AFT(3,:) = nanvar(MAT(LEARNT:size(Infos,1),:))./nanmean(MAT(LEARNT:size(Infos,1),:));
% % %         
% % %         
% % %         F = figure
% % %         % subplot(4,4,[6 7 10 11])
% % %         hold on;
% % %         % Align_code=11;
% % %         % Start = -750; End = 950;
% % %         X = Start:End; n=1;
% % %         
% % %         subplot(4,4,[1])
% % %         hold on;
% % %         Y = BEF(1,:); E = BEF(2,:);
% % %         colour1 = [205	112	84]/255;
% % %         errorline_n(X,Y,E,n,colour1,0.3)
% % %         xlim([Start+150 End-150])
% % %         YL = ylim;
% % %         ylim([-15 YL(2)])
% % %         
% % %         
% % %         subplot(4,4,[5])
% % %         hold on;
% % %         Y = DUR(1,:); E = DUR(2,:);
% % %         colour2 = [238	180	34]/255;
% % %         errorline_n(X,Y,E,n,colour2,0.3)
% % %         xlim([Start+150 End-150])
% % %         YL = ylim;
% % %         ylim([-15 YL(2)])
% % %         ylabel('Sp/s')
% % %         
% % %         subplot(4,4,[9])
% % %         hold on;
% % %         Y = AFT(1,:); E = AFT(2,:);
% % %         colour3 = [205	133	63]/255;
% % %         errorline_n(X,Y,E,n,colour3,0.3)
% % %         xlim([Start+150 End-150])
% % %         YL = ylim;
% % %         ylim([-15 YL(2)])
% % %         % xlabel('Time (ms)')
% % %         
% % %         subplot(4,4,[13])
% % %         hold on;
% % %         plot(X,BEF(2,:),'color',colour1,'linewidth',1.2)
% % %         plot(X,DUR(2,:),'color',colour2,'linewidth',1.2)
% % %         plot(X,AFT(2,:),'color',colour3,'linewidth',1.2)
% % %         xlim([Start+150 End-150])
% % %         % ylim([5 30])
% % %         % xlabel('Time (ms)')
% % %         
% % %         suptitle(strcat(NOME(1:8),'-',NOME(10:12)));
% % %         cd(Results_dir)
% % %         filename = strcat(NOME,'PSTH_CONDITIONS_',EPOCH_TYPE)
% % %         print(F, '-dpdf', filename, '-r600')
% % %         
% % %         
% % %         % note the time of max variance and the align code
% % %         if Align_code==4
% % %             [MAX,MPOS]=nanmax(DUR(2,find(X==Start+150):find(X==400)))
% % %             VAR_MAX(1) = X(MPOS)+150;
% % %             VAR_MAX(2) = Align_code
% % %         end
% % %         
% % %         if strcmp(EPOCH_TYPE,'P')
% % %             [MAX,MPOS]=nanmax(DUR(2,find(X==200):find(X==End-50)))
% % %             VAR_MAX(1) = X(MPOS)+find(X==200);
% % %             VAR_MAX(2) = Align_code
% % %         end
% % %         
% % %         if strcmp(EPOCH_TYPE,'M')
% % %             [MAX,MPOS]=nanmax(DUR(2,find(X==-300):find(X==End-550)))
% % %             VAR_MAX(1) = X(MPOS)+find(X==-300);
% % %             VAR_MAX(2) = Align_code
% % %         end
% % %         
% % %         
% % %         
% % %         
% % %         save(MERGE_file,'VAR_MAX','-append');
% % %         save(ALLCELLS_file,'VAR_MAX','-append');
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         %% Varience analysis across trials
% % %         
% % %         
% % %         %
% % %         %
% % %         % if Start>=500 & Align_code == 11
% % %         %     a = -1200;
% % %         %     b = 200;
% % %         %     End_RANDOM = round((b-a).*rand(1) + a);
% % %         %     Start_RANDOM = End_RANDOM-200;
% % %         % end
% % %         %
% % %         % if Start<=0 & Align_code == 4
% % %         %     a = 400;
% % %         %     b = 1200;
% % %         %     Start_RANDOM = round((b-a).*rand(1) + a);
% % %         %     End_RANDOM = Start_RANDOM+200;
% % %         % end
% % %         %
% % %         %
% % %         % if Start<=300 & Align_code == 4
% % %         %     a = 800;
% % %         %     b = 1200;
% % %         %     Start_RANDOM = round((b-a).*rand(1) + a);
% % %         %     End_RANDOM = Start_RANDOM+200;
% % %         % end
% % %         %
% % %         %
% % %         % if Start>=-300 & Align_code == 11
% % %         %     a = -1200;
% % %         %     b = -500;
% % %         %     End_RANDOM = round((b-a).*rand(1) + a);
% % %         %     Start_RANDOM = End_RANDOM-200;
% % %         % end
% % %         %
% % %         %
% % %         % Start = EPOCH_START;
% % %         % End = EPOCH_END;
% % %         %
% % %         
% % %         %%
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         Signal1 = mean(Matrix1,2);
% % %         
% % %         F = figure;
% % %         subplot(3,3,[1 2 3])
% % %         hold on;
% % %         plot(Signal1(:,1));
% % %         
% % %         for iiii=1:size(Infos,1)
% % %             if Infos(iiii,10)==2
% % %                 plot(iiii,Signal1(iiii)+5,'vk','markersize',2,'color','r')
% % %             end
% % %         end
% % %         xlim([1 size(Infos,1)])
% % %         xlabel('trial number','fontsize',8);
% % %         ylabel('mean Sp/s','fontsize',8);
% % %         plot([CHANGE CHANGE],ylim,'-k','linewidth',1.5,'color',[0.3 0.3 0.3])
% % %         plot([LEARNT LEARNT],ylim,'-k','linewidth',1.5,'color',[0.3 0.3 0.3])
% % %         
% % %         
% % %         ind_c = find(Infos(:,10)==1);
% % %         CORRECT_vals(1) = nanmean(Signal1(ind_c));
% % %         CORRECT_vals(2) = nanstd(Signal1(ind_c));
% % %         ind_w = find(Infos(:,10)==2);
% % %         WRONG_vals(1) = nanmean(Signal1(ind_w));
% % %         WRONG_vals(2) = nanstd(Signal1(ind_w));
% % %         
% % %         
% % %         subplot(3,3,4)
% % %         beeswarm_plot_n({Signal1(ind_c) Signal1(ind_w)},'color',[0 0.7 0.7; 1 0.5 0.5],'stats','none','size',2);
% % %         value = [1, 2];
% % %         set(gca, 'XTick', value);
% % %         set(gca, 'XTickLabel', {'Corr','Wrng'},'fontsize',10);
% % %         yh = ylabel('Sp/s');
% % %         title('All trials','fontsize',8);
% % %         
% % %         
% % %         ind_C = find(Infos(CHANGE:LEARNT,10)==1)+CHANGE-1;
% % %         ind_W = find(Infos(CHANGE:LEARNT,10)==2)+CHANGE-1;
% % %         subplot(3,3,5)
% % %         beeswarm_plot_n({Signal1(ind_C) Signal1(ind_W)},'color',[0 0.7 0.7; 1 0.5 0.5],'stats','none','size',2);
% % %         value = [1, 2];
% % %         set(gca, 'XTick', value);
% % %         set(gca, 'XTickLabel', {'Corr','Wrng'},'fontsize',10);
% % %         yh = ylabel('Sp/s');
% % %         title('Learning','fontsize',8);
% % %         
% % %         
% % %         ind_C_BEF = find(Infos(1:CHANGE-1,10)==1);
% % %         ind_C_DUR = find(Infos(CHANGE:LEARNT,10)==1)+CHANGE-1;
% % %         ind_C_AFT = find(Infos(LEARNT+1:size(Infos,1),10)==1)+LEARNT;
% % %         subplot(3,3,6)
% % %         hold on;
% % %         
% % %         x_plot = [1 2 3];
% % %         y_plot = [nanmean(Signal1(ind_C_BEF)) nanmean(Signal1(ind_C_DUR)) nanmean(Signal1(ind_C_AFT))];
% % %         e_plot = [nanstd(Signal1(ind_C_BEF)) nanstd(Signal1(ind_C_DUR)) nanstd(Signal1(ind_C_AFT))];
% % %         
% % %         bar(x_plot,y_plot,'facecolor',[0 0.7 0.7]);
% % %         e = errorbar(x_plot,y_plot,e_plot,'.');
% % %         e.Color = [0 0 0];
% % %         
% % %         p_BD = ranksum(Signal1(ind_C_BEF),Signal1(ind_C_DUR));
% % %         p_DA = ranksum(Signal1(ind_C_AFT),Signal1(ind_C_DUR));
% % %         
% % %         for pp=1:2
% % %             if pp==1 p=p_BD; end
% % %             if pp==2 p=p_DA; end
% % %             
% % %             if p<0.05 & p>0.01
% % %                 STAR='*';
% % %             elseif p<0.01 & p>0.001
% % %                 STAR='**';
% % %             elseif p<0.0001
% % %                 STAR='***';
% % %             else
% % %                 STAR='ns';
% % %             end
% % %             
% % %             if pp==1 STAR_B = STAR; end
% % %             if pp==2 STAR_A = STAR; end
% % %             
% % %         end
% % %         
% % %         YL = ylim;
% % %         text(1.2,YL(2)-10,STAR_B,'fontsize',10);
% % %         text(2.2,YL(2)-10,STAR_A,'fontsize',10);
% % %         
% % %         value = [1,2,3];
% % %         set(gca, 'XTick', value);
% % %         set(gca, 'XTickLabel', {'Bef','Dur','Aft'},'fontsize',10);
% % %         yh = ylabel('Sp/s');
% % %         title('Correct','fontsize',8);
% % %         
% % %         suptitle(strcat(NOME(1:8),'-',NOME(10:12),'-VARIANCE-',EPOCH_TYPE));
% % %         cd(Results_dir)
% % %         filename = strcat(NOME,'_SS_VARIANCE_Trial_by_Trail_',EPOCH_TYPE)
% % %         print(F, '-dpdf', filename, '-r600')
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         % Signal2 = mean(Matrix2,2);
% % %         W_bin = 0.1; S_bin=0.05;
% % %         [x_trial, SIG_trial1, per_corr] = LEARNING_CURVE_n(Signal1, CHANGE, Infos, W_bin, S_bin,'p');
% % %         % [x_trial, SIG_trial2, per_corr] = LEARNING_CURVE_n(Signal2, CHANGE, Infos, W_bin, S_bin,'p');
% % %         
% % %         SIG_trial1 = SIG_trial1(:,2).*SIG_trial1(:,2);
% % %         % SIG_trial2 = SIG_trial2(:,2).*SIG_trial2(:,2);
% % %         
% % %         
% % %         
% % %         
% % %         clear VAR_ACTUAL VAR_RANDOM
% % %         ind = find(x_trial<CHANGE);
% % %         VAR_ACTUAL(1,1)=nanmean([SIG_trial1(ind(length(ind))) SIG_trial1(ind(length(ind)-1))]);
% % %         % VAR_ACTUAL(1,1)=nanmean([SIG_trial1(1) SIG_trial1(2)]);
% % %         % VAR_RANDOM(1,1)=nanmean([SIG_trial2(ind(length(ind))) SIG_trial2(ind(length(ind)-1))]);
% % %         
% % %         VAR_ACTUAL(2,1)=nanmean(SIG_trial1(find(x_trial>=CHANGE,3)));
% % %         % VAR_RANDOM(2,1)=nanmean(SIG_trial2(find(x_trial>=CHANGE,3)));
% % %         
% % %         VAR_ACTUAL(3,1)=nanmean([SIG_trial1(length(x_trial)) SIG_trial1(length(x_trial)-1)]);
% % %         % VAR_RANDOM(3,1)=nanmean([SIG_trial2(length(x_trial)) SIG_trial2(length(x_trial)-1)]);
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         F = figure
% % %         WID = 12;
% % %         LEN = 3;
% % %         FS=10;
% % %         
% % %         ax1 = subplot(WID,LEN,[1 2 4 5 7 8 13 14 16 17 19 20 22 23]);
% % %         hold on;
% % %         
% % %         plot(x_trial,per_corr,'-','color',[0.6 0.6 0.6],'linewidth',1)
% % %         plot(x_trial,per_corr,'Marker','O','MarkerSize',4,'Markerfacecolor',[0.5 0.5 0.5],'color',[0.5 0.5 0.5])
% % %         
% % %         % [H,H11,H21] = plotyy(x_trial,NaN(length(x_trial),1),x_trial,SIG_trial2(:,1));
% % %         % % [h,H1,H2] = plotyy(x_trial,NaN(length(x_trial),1),x_trial,SIG_trial2(:,1));
% % %         % % set(H2,'color',[169 130 181]/255,'linewidth',1,'Marker','^','MarkerSize',4,'Markerfacecolor',[169 130 181]/255,'color',[169 130 181]/255)
% % %         
% % %         
% % %         % [H,H1,H2] = plotyy(x_trial,per_corr,x_trial,SIG_trial1(:,1));
% % %         [H,H3,H4] = plotyy(x_trial,per_corr,x_trial,SIG_trial1(:,1));
% % %         set(H3,'linewidth',1,'Marker','O','MarkerSize',4,'Markerfacecolor',[0.5 0.5 0.5],'color',[0.5 0.5 0.5])
% % %         set(H4,'linewidth',1,'Marker','O','MarkerSize',4,'Markerfacecolor',[169 130 181]/255,'color',[169 130 181]/255)
% % %         set(H,{'ycolor'},{[0.6 0.6 0.6];[169 130 181]/255})
% % %         
% % %         plot([CHANGE CHANGE],ylim,'--','color',[0.6 0.6 0.6])
% % %         box off;
% % %         
% % %         
% % %         set(H(2),'YLim',[0 nanmax(SIG_trial1(:,1))])
% % %         % set(h(2),'YLim',[0 nanmax(SIG_trial1(:,1))])
% % %         
% % %         
% % %         
% % %         
% % %         % YL = ylim;
% % %         % plot(CHANGE,YL(1),'^k','MarkerFaceColor','K')
% % %         xlabel('trial number');
% % %         ylabel(H(1),'% Corr')
% % %         ylabel(H(2),'Variance (sp/s^2)')
% % %         
% % %         ax3 = subplot(WID,LEN,[6 9 12 15]);
% % %         hold on;
% % %         % plot(VAR_RANDOM,'color',[169 130 181]/255)
% % %         % plot(VAR_RANDOM,'^','color',[169 130 181]/255)
% % %         plot(VAR_ACTUAL,'color',[169 130 181]/255)
% % %         plot(VAR_ACTUAL,'o','color',[169 130 181]/255)
% % %         box off;
% % %         
% % %         
% % %         
% % %         [coeff1,p1] = corr(SIG_trial1(:,1),per_corr)
% % %         % [coeff2,p2] = corr(SIG_trial2(:,1),per_corr)
% % %         
% % %         
% % %         
% % %         ax4 = subplot(WID,LEN,[18 21 24 27 30]);
% % %         hold on;
% % %         axis off;
% % %         text(0,0.7,strcat('Actual = ',num2str(EPOCH_START),':',num2str(EPOCH_END),' @',ALIGN),'fontsize',10);
% % %         % text(0,0.6,strcat('Random = ',num2str(Start_RANDOM),':',num2str(End_RANDOM),' @',ALIGN),'fontsize',10);
% % %         text(0,0.4,strcat('RHO = ',num2str(coeff1)),'fontsize',10);
% % %         text(0,0.3,strcat('p = ',num2str(p1)),'fontsize',10);
% % %         % text(0,0.1,strcat('RHO RANDOM =',num2str(coeff2)),'fontsize',10);
% % %         % text(0,0,strcat('p RANDOM = ',num2str(p2)),'fontsize',10);
% % %         
% % %         
% % %         subplot(WID,LEN,[31 32])
% % %         imagesc(per_corr(:,1)')
% % %         colormap(jet)
% % %         axis off;
% % %         
% % %         clear TEMP_SIG
% % %         TEMP_SIG(:,1) = SIG_trial1(:,1);
% % %         % TEMP_SIG(:,2) = SIG_trial2(:,1);
% % %         TEMP_SIG = mat2gray(TEMP_SIG);
% % %         
% % %         subplot(WID,LEN,[34 35])
% % %         imagesc(TEMP_SIG')
% % %         colormap(jet)
% % %         axis off;
% % %         
% % %         
% % %         suptitle(strcat(NOME(1:8),'-',NOME(10:12),'-VARIANCE-',EPOCH_TYPE));
% % %         
% % %         cd(Results_dir)
% % %         filename = strcat(NOME,'_SS_VARIANCE_',EPOCH_TYPE)
% % %         print(F, '-dpdf', filename, '-r600')
% % %         
% % %         
% % %         SIG_VAR_XVAL = x_trial;
% % %         SIG_VAR_PCOR = per_corr;
% % %         SIG_VAR_TRUE = SIG_trial1;
% % %         % SIG_VAR_RAND = SIG_trial2;
% % %         
% % %         RHO_VAR_TRUE = coeff1;
% % %         p_VAR_TRUE = p1;
% % %         % RHO_VAR_RAND = coeff2;
% % %         % p_VAR_RAND = p2;
% % %         
% % %         SIG_VAR_XVAL = x_trial;
% % %         SIG_VAR_PCOR = per_corr;
% % %         SIG_VAR_TRUE = SIG_trial1;
% % %         % SIG_VAR_RAND = SIG_trial2;
% % %         
% % %         
% % %         
% % %         %
% % %         % cd(PathName1)
% % %         % save(filenm,'RHO_VAR_TRUE','p_VAR_TRUE','RHO_VAR_RAND','p_VAR_RAND','VAR_ACTUAL','VAR_RANDOM','-append');
% % %         % save(DATA_file,'RHO_VAR_TRUE','p_VAR_TRUE','RHO_VAR_RAND','p_VAR_RAND','VAR_ACTUAL','VAR_RANDOM','-append');
% % %         % save(POP_file,'RHO_VAR_TRUE','p_VAR_TRUE','RHO_VAR_RAND','p_VAR_RAND','VAR_ACTUAL','VAR_RANDOM','-append');
% % %         
% % %         
% % %         
% % %         save(MERGE_file,'RHO_VAR_TRUE','p_VAR_TRUE','VAR_ACTUAL','-append');
% % %         save(ALLCELLS_file,'RHO_VAR_TRUE','p_VAR_TRUE','VAR_ACTUAL','-append');
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         
% % %         VARIANCE.RHO(index,1) = RHO_VAR_TRUE;
% % %         VARIANCE.P(index,1) = p_VAR_TRUE;
% % %         VARIANCE.values{index,1} = VAR_ACTUAL;
% % %         VARIANCE.max(index,1) = VAR_MAX(1);
% % %         
% % %         
% % %     end
% % % end
% % % 


%% STATE CHANGE ANALYSIS PROPER


Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;


F = figure();

subplot(3,3,1)
hold on;
IND = 1:CHANGE;
P_MAIN_T = PSTHe_n(Spikes.S(IND),Infos(IND,4),Start_T,End_T,30,[0.5 0.5 0.5],1,1);
temp_time = Start_T:End_T;
temp_ind = find(StartLIM_T<=temp_time & temp_time<=EndLIM_T);
errorline_n(StartLIM_T:EndLIM_T,P_MAIN_T(1,temp_ind),P_MAIN_T(3,temp_ind),1,[0.5 0.5 0.5],0.3,0,1);

IND = LEARNT:size(Infos,1);
P_LEARNT_T = PSTHe_n(Spikes.S(IND),Infos(IND,4),Start_T,End_T,30,[0.5 0.5 0.5],1,1);
temp_time = Start_T:End_T;
temp_ind = find(StartLIM_T<=temp_time & temp_time<=EndLIM_T);
errorline_n(StartLIM_T:EndLIM_T,P_LEARNT_T(1,temp_ind),P_LEARNT_T(3,temp_ind),1,[0 0 0],0.3,0,1);

xlim([StartLIM_T EndLIM_T]);

MAIN_NUM = CHANGE-1;
AFTER_NUM = size(Infos,1)-LEARNT;

%%%% stats
for st=1:length(temp_ind)
    p(1,st)=ttest_n(P_MAIN_T(1,st),P_LEARNT_T(1,st),P_MAIN_T(2,st),P_LEARNT_T(2,st),MAIN_NUM,AFTER_NUM);
end
IND_p{1,1} = find(p(1,:)<0.05);
YLIM=ylim;
plot_time=Start_T:End_T;
plot(plot_time(IND_p{1,1}),[YLIM(2)*ones(length(IND_p{1,1}),1)],'.','color',[212 175 55]/255);
xlim([StartLIM_T EndLIM_T]);


subplot(3,3,2)
hold on;
IND = 1:CHANGE;
P_MAIN_M = PSTHe_n(Spikes.S(IND),Infos(IND,11),Start_M,End_M,30,[0.5 0.5 0.5],1,1);
temp_time = Start_M:End_M;
temp_ind = find(StartLIM_M<=temp_time & temp_time<=EndLIM_M);
errorline_n(StartLIM_M:EndLIM_M,P_MAIN_M(1,temp_ind),P_MAIN_M(3,temp_ind),1,[0.5 0.5 0.5],0.3,0,1);

IND = LEARNT:size(Infos,1);
P_LEARNT_M = PSTHe_n(Spikes.S(IND),Infos(IND,11),Start_M,End_M,30,[0.5 0.5 0.5],1,1);
temp_time = Start_M:End_M;
temp_ind = find(StartLIM_M<=temp_time & temp_time<=EndLIM_M);
errorline_n(StartLIM_M:EndLIM_M,P_LEARNT_M(1,temp_ind),P_LEARNT_M(3,temp_ind),1,[0 0 0],0.3,0,1);

xlim([StartLIM_M EndLIM_M]);


MAIN_NUM = CHANGE-1;
AFTER_NUM = size(Infos,1)-LEARNT;

%%%% stats
for st=1:length(temp_ind)
    p(2,st)=ttest_n(P_MAIN_M(1,st),P_LEARNT_M(1,st),P_MAIN_M(2,st),P_LEARNT_M(2,st),MAIN_NUM,AFTER_NUM);
end
IND_p{2,1} = find(p(2,:)<0.05);
YLIM=ylim;
plot_time=Start_M:End_M;
plot(plot_time(IND_p{2,1}),[YLIM(2)*ones(length(IND_p{2,1}),1)],'.','color',[212 175 55]/255);

xlim([StartLIM_M EndLIM_M]);





bin=round(0.20*size(Infos,1));
% beginColor = [193 188 118]/255;
% endColor = [114 176 149]/255;

beginColor = [11 84 42]/255;
endColor = [138 226 194]/255;


numSteps = length(1:bin:size(Infos,1)-bin);
cMap = makeColorMap(beginColor, endColor, numSteps);


subplot(3,3,4)
hold on;
COUNT=0;
for i=1:bin:size(Infos,1)-bin
    IND = i:i+bin;
    COUNT = COUNT+1;
    P_MAIN_T = PSTH_n(Spikes.S(IND),Infos(IND,4),Start_T,End_T,30,cMap(COUNT,:),0.7,0);
end

xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])


subplot(3,3,5)
hold on;
COUNT=0;
for i=1:bin:size(Infos,1)-bin
    IND = i:i+bin;
    COUNT = COUNT+1;
    P_MAIN_M = PSTH_n(Spikes.S(IND),Infos(IND,11),Start_M,End_M,30,cMap(COUNT,:),0.7,0);
end

xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])



subplot(3,3,6)
surf(zeros(2));
colormap(cMap)
colorbar
axis off;



suptitle(strcat('STATE-',NOME(1:8),'-',NOME(10:12)));


cd(Results_dir)
filename = strcat(NOME,'_STATE_CHANGE');
print(F, '-dpdf', filename, '-r400')






% % % 
% % % 
% % % for index=1:4
% % %     STATE.FR(index,2) = NaN;
% % % end
% % % 
% % % for index=1:2
% % %     if DELTA.EpochFlag(index)==1
% % %         Time = Start_T:End_T;
% % %         INDS = find(DELTA.START(index)+75<=Time & Time<=DELTA.END(index)-75);
% % %         STATE.FR(index,1) = nanmean(P_MAIN_T(1,INDS));
% % %         STATE.FR(index,2) = nanmean(P_LEARNT_T(1,INDS));
% % %     end
% % % end
% % % 
% % % for index=3:4
% % %     if DELTA.EpochFlag(index)==1
% % %         Time = Start_M:End_M;
% % %         INDS = find(DELTA.START(index)+75<=Time & Time<=DELTA.END(index)-75);
% % %         STATE.FR(index,1) = nanmean(P_MAIN_M(1,INDS));
% % %         STATE.FR(index,2) = nanmean(P_LEARNT_M(1,INDS));
% % %     end
% % % end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONS ANALYSIS %%

COLOUR = [218 165 32]/255;
% Trial by chance
TBC_W = find(Infos(CHANGE+1:size(Infos,1),10)==2,1)+CHANGE;
TBC_C = find(Infos(CHANGE+1:size(Infos,1),10)==1,1)+CHANGE;

SIGMA = 50;

F = figure();
subplot(4,4,6)
PSTH_n(Spikes.S(CHANGE-1),Infos(CHANGE-1,4),Start_T,End_T,SIGMA,[0.2 0.2 0.2],0.7,0);
PSTH_n(Spikes.S(CHANGE),Infos(CHANGE,4),Start_T,End_T,SIGMA,[0.6 0.6 0.6],0.7,0);
% PSTH_n(Spikes.S(TBC_W),Infos(TBC_W,4),Start_T,End_T,SIGMA,'r',0.7,0);
% PSTH_n(Spikes.S(TBC_C),Infos(TBC_C,4),Start_T,End_T,SIGMA,'b',0.7,0);
PSTH_n(Spikes.S(CHANGE+1),Infos(CHANGE+1,4),Start_T,End_T,SIGMA,COLOUR,0.7,0);
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])
YLIM1 = ylim;


subplot(4,4,7)
PSTH_n(Spikes.S(CHANGE-1),Infos(CHANGE-1,11),Start_M,End_M,SIGMA,[0.2 0.2 0.2],0.7,0);
PSTH_n(Spikes.S(CHANGE),Infos(CHANGE,11),Start_M,End_M,SIGMA,[0.6 0.6 0.6],0.7,0);
% PSTH_n(Spikes.S(TBC_W),Infos(TBC_W,11),Start_M,End_M,SIGMA,'r',0.7,0);
% PSTH_n(Spikes.S(TBC_C),Infos(TBC_C,11),Start_M,End_M,SIGMA,'b',0.7,0);
PSTH_n(Spikes.S(CHANGE+1),Infos(CHANGE+1,11),Start_M,End_M,SIGMA,COLOUR,0.7,0);
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])
YLIM2 = ylim;

MIN = nanmin(YLIM1(1),YLIM2(1));
MAX= nanmax(YLIM1(2),YLIM2(2));

subplot(4,4,6)
ylim([MIN MAX]);
subplot(4,4,7)
ylim([MIN MAX]);


% FIND THE FIRST TRIAL AFTER CHANGE -------------------------

if Infos(CHANGE-1,10)==1 COLOUR1 = 'b'; end
if Infos(CHANGE-1,10)==2 COLOUR1 = 'r'; end

if Infos(CHANGE,10)==1 COLOUR2 = 'b'; end
if Infos(CHANGE,10)==2 COLOUR2 = 'r'; end

if Infos(CHANGE+1,10)==1 COLOUR3 = 'b'; end
if Infos(CHANGE+1,10)==2 COLOUR3 = 'r'; end

if Infos(CHANGE+2,10)==1 COLOUR4 = 'b'; end
if Infos(CHANGE+2,10)==2 COLOUR4 = 'r'; end

if Infos(CHANGE+3,10)==1 COLOUR5 = 'b'; end
if Infos(CHANGE+3,10)==2 COLOUR5 = 'r'; end

subplot(4,4,10)
% for i=CHANGE-5:CHANGE-1
%     PSTH_ONE_n(Spikes.S{i},Infos(i,4),Start_T,End_T,SIGMA,[0.7 0.7 0.7],0.7,0); 
% end
i=CHANGE-6:CHANGE-2;
PSTHe_n(Spikes.S(i),Infos(i,4),Start_T,End_T,SIGMA,[0.7 0.7 0.7],0.7,0); 

set(gca,'FontSize',2,'LineWidth',0.7)
PSTH_n(Spikes.S(CHANGE-1),Infos(CHANGE-1,4),Start_T,End_T,SIGMA,COLOUR1,0.2,0);
PSTH_n(Spikes.S(CHANGE),Infos(CHANGE,4),Start_T,End_T,SIGMA,COLOUR2,1,0);
PSTH_n(Spikes.S(CHANGE+1),Infos(CHANGE+1,4),Start_T,End_T,SIGMA,COLOUR3,2,0);
PSTH_n(Spikes.S(CHANGE+2),Infos(CHANGE+2,4),Start_T,End_T,SIGMA,COLOUR4,3,0);
PSTH_n(Spikes.S(CHANGE+3),Infos(CHANGE+3,4),Start_T,End_T,SIGMA,COLOUR5,4,0);
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])



subplot(4,4,11)
% for i=CHANGE-5:CHANGE-1
%     PSTH_ONE_n(Spikes.S{i},Infos(i,11),Start_M,End_M,SIGMA,[0.7 0.7 0.7],0.7,0); 
% end

i=CHANGE-6:CHANGE-2;
PSTHe_n(Spikes.S(i),Infos(i,11),Start_M,End_M,SIGMA,[0.7 0.7 0.7],0.7,0); 

set(gca,'FontSize',2,'LineWidth',0.7)
PSTH_n(Spikes.S(CHANGE-1),Infos(CHANGE-1,11),Start_M,End_M,SIGMA,COLOUR1,0.2,0);
PSTH_n(Spikes.S(CHANGE),Infos(CHANGE,11),Start_M,End_M,SIGMA,COLOUR2,1,0);
PSTH_n(Spikes.S(CHANGE+1),Infos(CHANGE+1,11),Start_M,End_M,SIGMA,COLOUR3,2,0);
PSTH_n(Spikes.S(CHANGE+2),Infos(CHANGE+2,11),Start_M,End_M,SIGMA,COLOUR4,3,0);
PSTH_n(Spikes.S(CHANGE+3),Infos(CHANGE+3,11),Start_M,End_M,SIGMA,COLOUR5,4,0);
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])




suptitle(strcat('CONS-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = strcat(NOME,'_STATE_CHANGE_CONS');
print(F, '-dpdf', filename, '-r400')















%%%%%%%%%%%%% all corr and wrong trials during learning

% ENDTRIAL = LEARNT;
ENDTRIAL = CHANGE+20;

ALL_W = find(Infos(CHANGE-1:ENDTRIAL,10)==2)+CHANGE-2;
ALL_C = find(Infos(CHANGE-1:ENDTRIAL,10)==1)+CHANGE-2;



beginColor = [76 5 5]/255;
endColor = [255 180 180]/255;
numSteps = length(ALL_W);
W_COLOUR = makeColorMap(beginColor, endColor, numSteps);


beginColor = [5 54 5]/255;
endColor = [193 247 193]/255;
numSteps = length(ALL_C);
C_COLOUR = makeColorMap(beginColor, endColor, numSteps);


F = figure();

subplot(4,4,6)
for i=1:length(ALL_W)
    PSTH_ONE_n(Spikes.S{ALL_W(i)},Infos(ALL_W(i),4),Start_T,End_T,SIGMA,W_COLOUR(i,:),0.7,0); 
end
set(gca,'FontSize',2,'LineWidth',0.7)
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])


subplot(4,4,7)
for i=1:length(ALL_W)
    PSTH_ONE_n(Spikes.S{ALL_W(i)},Infos(ALL_W(i),11),Start_M,End_M,SIGMA,W_COLOUR(i,:),0.7,0); 
end
set(gca,'FontSize',2,'LineWidth',0.7)
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])




subplot(4,4,10)
for i=1:length(ALL_C)
    PSTH_ONE_n(Spikes.S{ALL_C(i)},Infos(ALL_C(i),4),Start_T,End_T,SIGMA,C_COLOUR(i,:),0.7,0); 
end
set(gca,'FontSize',2,'LineWidth',0.7)
xlim([StartLIM_T EndLIM_T]);
xlabel([])
ylabel([])


subplot(4,4,11)
for i=1:length(ALL_C)
    PSTH_ONE_n(Spikes.S{ALL_C(i)},Infos(ALL_C(i),11),Start_M,End_M,SIGMA,C_COLOUR(i,:),0.7,0); 
end
set(gca,'FontSize',2,'LineWidth',0.7)
xlim([StartLIM_M EndLIM_M]);
xlabel([])
ylabel([])


suptitle(strcat('CONS-CW-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = strcat(NOME,'_STATE_CHANGE_CONS_CW');
print(F, '-dpdf', filename, '-r400')



% % % save(POP_file,'STATE','VARIANCE','-append');
% % % save(MERGE_file,'STATE','VARIANCE','-append');
% % % save(ALLCELLS_file,'STATE','VARIANCE','-append');








STATE_FIGURE_ALL






% end
