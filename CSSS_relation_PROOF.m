

%% CSSS_relation_PROOF
% Created by Naveen on 8/4/2019 at JLG


C_SPK = Spikes.C;
S_SPK = Spikes.S; 
Start=-50; End = 50;

try
    C_SPKshape = SPK_align_dtw(cell2mat(Shapes.C),1);
    S_SPKshape = SPK_align_dtw(cell2mat(Shapes.S),1);
catch
    C_SPKshape = cell2mat(Shapes.C);
    S_SPKshape = cell2mat(Shapes.S);
end

F = figure();


%% waveforms ---------------------
s1 = subplot(4,4,1);
hold on;
errorline_n(1:size(S_SPKshape,2),nanmean(S_SPKshape),nanstd(S_SPKshape)/sqrt(size(S_SPKshape,1)),1,[0 0 0],0.3,0,0.6);
xlim([1 size(S_SPKshape,2)]);
% ylim([-1 1]);
YLIM = ylim;
plot([0 50],[YLIM(2) YLIM(2)],'-k');
text(0, -0.8,'1 ms','fontsize',6.7);
axis off;

s2 = subplot(4,4,5);
hold on;
errorline_n(1:size(C_SPKshape,2),nanmean(C_SPKshape),nanstd(C_SPKshape)/sqrt(size(C_SPKshape,1)),1,[0 0 0],0.3,0,0.6);
xlim([1 size(C_SPKshape,2)]);
% ylim([-1 1]);
YLIM = ylim;
plot([0 50],[YLIM(2) YLIM(2)],'-k');
text(0, -0.8,'1 ms','fontsize',6.7);
axis off;

linkaxes([s1 s2],'xy');


%% SS_ISI -------------------------
s1 = subplot(4,4,2);
hold on;

clear temp TEMP TEMP_new
cellflag=0;
parfor i=1:size(S_SPK,1)
    Signal_S = S_SPK{i,1};
    cell_len = length(S_SPK{i,1});
    if cell_len>=2
        temp{i,1} = diff(Signal_S);
        cellflag=cellflag+1;
    end
end

if cellflag>0
    
    TEMP = sort(cell2mat(temp));
    TEMP_new = TEMP;
    
    BIN_W = 1;
    
%     hh = histcounts(TEMP_new,'BinWidth',BIN_W)/length(TEMP_new);
%     plot(hh);
    
    h = histogram(TEMP_new,'BinWidth',BIN_W,'DisplayStyle','stairs');
    
    YLIM = round(nanmax(h.BinCounts)/length(TEMP_new),2);

    yticks([0 nanmax(h.BinCounts)])
    yticklabels({'0',num2str(YLIM)})
    
    if ~isempty(TEMP_new)
%         ylim([0 YLIM*length(TEMP_new)]);
%         YLIMM(1) = YLIM*length(TEMP_new);
        xlim([0 50]);
    end
    set(gca,'FontSize',8,'LineWidth',0.2)
    ylabel([]);
end

    
%% CS_ISI -------------------------
s2 = subplot(4,4,6);
hold on;

clear temp TEMP TEMP_new
cellflag=0;
parfor i=1:size(C_SPK,1)
    Signal_C = C_SPK{i,1};
    cell_len = length(C_SPK{i,1});
    if cell_len>=2
        temp{i,1} = diff(Signal_C);
        cellflag=cellflag+1;
    end
end

if cellflag>0
    TEMP = sort(cell2mat(temp));
    TEMP_new = TEMP;
    
    BIN_W = 100;
    h = histogram(TEMP_new,'BinWidth',BIN_W,'DisplayStyle','stairs');
%     h.FaceColor = [0 0 0];
   
    YLIM = round(nanmax(h.BinCounts)/length(TEMP_new),2);
    
    yticks([0 nanmax(h.BinCounts)])
    yticklabels({'0',num2str(YLIM)})
    
    if ~isempty(TEMP_new)
%         ylim([0 YLIM*length(TEMP_new)]);
%         YLIMM(2) = YLIM*length(TEMP_new);
        xlim([0 5000]);
    end
    set(gca,'FontSize',8,'LineWidth',0.2)
    ylabel([]);
end   
    
%    
% linkaxes([s1,s2],'y')
% nanmax(YLIMM)


subplot(4,4,3)
hold on;

% P(S(t)|S(0)) -------------------------
S_S = CER_SS_PSTH_ALIGNED(S_SPK,Start-50,End+50,1,[0 0 1]);

% P(S(t)|C(0)) -------------------------
C_S = CER_CS_PSTH_ALIGNED(S_SPK,C_SPK,Start-50,End+50,1,[1 0 0]);

% P(S(t)|C(x)) -------------------------
r = randi([1,length(S_SPK)],length(S_SPK),1);
S_SPK = S_SPK(r); % shuffle all the concatenated S spikes from 3 random cells
S_SPK = S_SPK(1:length(C_SPK)); % Take only that many trials as the main cell
CER_CS_PSTH_ALIGNED(S_SPK,C_SPK,Start,End,1,[0.5 0.5 0.5])

ylabel([]);
set(gca,'fontsize',8)
box off;
xlim([Start End])




filename = 'SS_CS_RELATION';
cd(Results_dir)
print(F, '-dpdf', filename, '-r400')










save(POP_file,'S_S','C_S','-append');







 
% end