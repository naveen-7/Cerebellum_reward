

F = figure
subplot (3,3,5);

clear temp TEMP TEMP_new
cellflag=0;
for i=1:size(Spikes.S,1)
    Signal_S = Spikes.S{i,1};
    Time_F = Infos(i,3);   % FP_on
    Time_S = Infos(i,4);   % Stimulus
    
    Length = Time_S-Time_F;
    
    Sig_S = Signal_S(find(Time_F<Signal_S & Signal_S<Time_S));
    Sig_S = Sig_S-Time_F;
    cell_len = length(Spikes.S{i,1});
    if cell_len>=2
        temp{i,1} = diff(Sig_S);
        cellflag=cellflag+1;
    end
end



if cellflag>0
    
    TEMP = sort(cell2mat(temp));
    
    %     temp_avg = nanmean(TEMP);
    %     temp_std = nanstd(TEMP);
    %     TEMP_new = TEMP(1:find(TEMP>=temp_avg+temp_std,1));
    
    TEMP_new = TEMP;
    
    BIN_W = 7;
    h = histogram(TEMP_new,'BinWidth',BIN_W);
    %     h = histogram(TEMP_new);
    h.FaceColor = [0.1255    0.6980    0.6667];
    
    xlabel('ISI (in ms)','FontSize',7);
    if ~isempty(TEMP_new)
        xlim([min(TEMP_new) max(TEMP_new)]);
        %         ylim([min(nb) max(nb)]);
    end
    
    title('ISI','FontSize',8);
    set(gca,'FontSize',7,'LineWidth',0.3)
    ylabel([]);
    box off;
    
    MED_ISI = nanmedian(TEMP_new);
    YLIM = ylim;
    text(max(TEMP_new)*0.6,YLIM(2)*0.7,strcat('ISI = ',num2str(MED_ISI)),'fontsize',7);
    
    
    cd(Results_dir)
    filename = strcat(NOME,'_isi');
    print(F, '-dpdf', filename, '-r400')
    
    
    
    
    
    save(MERGE_file,'MED_ISI','-append');
    save(ALLCELLS_file,'MED_ISI','-append');
    save(POP_file,'MED_ISI','-append');
    
    
end






