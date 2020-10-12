

%% function CER_analysis_03c_LINEAR_rsm
% Memory analysis
% written by naveen at JLG on 9/10/19

NUM = 19;
SIG = 40;

%%
% 
WRG_COLOUR = [255 128 128]/255;
COR_COLOUR = [0 179 179]/255;
OT_COLOUR = [81 81 81]/255;



COLOR = [168 096 143; ...   % WC; mismatch
         235 052 171; ...   % WC; match
         100 098 158; ...   % CW; mismatch
         000 000 255; ...   % CW; match
         097 140 102; ...   % CC; mismatch
         000 255 000; ...   % CC; match
         171 097 097; ...   % WW; mismatch
         255 000 000; ...   % WW; match
    ]/255;



Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;



%%%%%% CODE ALL THE SYMBOLS 
SYM(:,1) = Infos(:,5);  % symbol
SYM(find(SYM==888))=NaN;
SYM(find(SYM==999))=NaN;
UN = unique(SYM(CHANGE+1:end));
UN(:,2)=1:length(UN); 
for ii=1:length(UN)
    SYM(find(SYM==UN(ii,1)))=UN(ii,2);
end


%%%% RSM structure
%%% 1 Rew
%%% 2 Sym
%%% 3 RS
%%% 4 Mov
%%% 5 RM


clear RSM
RSM = zeros(length(Infos),5);
%%%% SEPERATE BASED ON REWARD
for i=2:length(Infos)
    
    %%%% REWARD
    if (Infos(i,10)==1) && (Infos(i-1,10)==2)  RSM(i,1)=1;  end   % Wrong -> Correct
    if (Infos(i,10)==2) && (Infos(i-1,10)==1)  RSM(i,1)=2;  end   % Correct -> Wrong 
    if (Infos(i,10)==1) && (Infos(i-1,10)==1)  RSM(i,1)=3;  end   % Correct -> Correct
    if (Infos(i,10)==2) && (Infos(i-1,10)==2)  RSM(i,1)=4;  end   % Wrong -> Wrong
    
    %%%% SYMBOL
    if (SYM(i,1)==1) && (SYM(i-1,1)==2)  RSM(i,2)=1;  end % 1 -> 2    
    if (SYM(i,1)==2) && (SYM(i-1,1)==1)  RSM(i,2)=1;  end % 2 -> 1    
    if (SYM(i,1)==1) && (SYM(i-1,1)==1)  RSM(i,2)=2;  end % 1 -> 1
    if (SYM(i,1)==2) && (SYM(i-1,1)==2)  RSM(i,2)=2;  end % 2 -> 2
    
    RSM(i,3) = RSM(i,1)*RSM(i,2);
    if (RSM(i,1)==2 && RSM(i,2)==1) RSM(i,3)=3; end;
    if (RSM(i,1)==3 && RSM(i,2)==1) RSM(i,3)=5; end;
    if (RSM(i,1)==4 && RSM(i,2)==1) RSM(i,3)=7; end;
    
    
    %%%% MOVEMENT
    if (Infos(i,9)==1) && (Infos(i-1,9)==2)  RSM(i,4)=1;  end % 1 -> 2    
    if (Infos(i,9)==2) && (Infos(i-1,9)==1)  RSM(i,4)=1;  end % 2 -> 1    
    if (Infos(i,9)==1) && (Infos(i-1,9)==1)  RSM(i,4)=2;  end % 1 -> 1
    if (Infos(i,9)==2) && (Infos(i-1,9)==2)  RSM(i,4)=2;  end % 2 -> 2
    
    RSM(i,5) = RSM(i,1)*RSM(i,4);
    if (RSM(i,1)==2 && RSM(i,4)==1) RSM(i,5)=3; end;
    if (RSM(i,1)==3 && RSM(i,4)==1) RSM(i,5)=5; end;
    if (RSM(i,1)==4 && RSM(i,4)==1) RSM(i,5)=7; end;
    
end

RSM(RSM==0)=NaN;







%%%%%%%%%%%%% SYMBOL %%%%%%%%%%%%%%%%%%%%%%%%%%%
LEN = zeros(8,1);
FF = figure();
subplot(2,2,1)
hold on;
for i=1:8
    clear IND;
    IND = find(RSM(CHANGE:CHANGE+NUM,3)==i)+CHANGE-1;
    Psym_T(i,:) = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,COLOR(i,:),1,0); % Aligned to target
    LEN(i,1) = length(IND);
end

subplot(2,2,2)
hold on;
for i=1:8
    clear IND;
    IND = find(RSM(CHANGE:CHANGE+NUM,3)==i)+CHANGE-1;
    Psym_M(i,:) =PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COLOR(i,:),1,0); % Aligned to target
end

subplot(2,2,3)
box off; axis off;
text(0,0.1,strcat('Wrong->Correct, nomatch =',num2str(LEN(1))),'color',COLOR(1,:),'fontsize',8);
text(0,0.2,strcat('Wrong->Correct, MATCH =',num2str(LEN(2))),'color',COLOR(2,:),'fontsize',8);
text(0,0.3,strcat('Correct->Wrong, nomatch =',num2str(LEN(3))),'color',COLOR(3,:),'fontsize',8);
text(0,0.4,strcat('Correct->Wrong, MATCH =',num2str(LEN(4))),'color',COLOR(4,:),'fontsize',8);
text(0,0.5,strcat('Correct->Correct, nomatch =',num2str(LEN(5))),'color',COLOR(5,:),'fontsize',8);
text(0,0.6,strcat('Correct->Correct, MATCH =',num2str(LEN(6))),'color',COLOR(6,:),'fontsize',8);
text(0,0.7,strcat('Wrong->Wrong, nomatch =',num2str(LEN(7))),'color',COLOR(7,:),'fontsize',8);
text(0,0.8,strcat('Wrong->Wrong, MATCH =',num2str(LEN(8))),'color',COLOR(8,:),'fontsize',8);


suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12),'-sym'));
cd(Results_dir)
filename = 'PSTH_Simple spike_OCTO_SYM';
print(FF, '-dpdf', filename, '-r400')






%%%%%%%%%%%%% MOVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%
LEN = zeros(8,1);
FF = figure();
subplot(2,2,1)
hold on;
for i=1:8
    clear IND;
    IND = find(RSM(CHANGE:CHANGE+NUM,5)==i)+CHANGE-1;
    Pmov_T(i,:) = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,COLOR(i,:),1,0); % Aligned to target
    LEN(i,1) = length(IND);
end

subplot(2,2,2)
hold on;
for i=1:8
    clear IND;
    IND = find(RSM(CHANGE:CHANGE+NUM,5)==i)+CHANGE-1;
    Pmov_M(i,:) =PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COLOR(i,:),1,0); % Aligned to target
end

subplot(2,2,3)
box off; axis off;
text(0,0.1,strcat('Wrong->Correct, nomatch =',num2str(LEN(1))),'color',COLOR(1,:),'fontsize',8);
text(0,0.2,strcat('Wrong->Correct, MATCH =',num2str(LEN(2))),'color',COLOR(2,:),'fontsize',8);
text(0,0.3,strcat('Correct->Wrong, nomatch =',num2str(LEN(3))),'color',COLOR(3,:),'fontsize',8);
text(0,0.4,strcat('Correct->Wrong, MATCH =',num2str(LEN(4))),'color',COLOR(4,:),'fontsize',8);
text(0,0.5,strcat('Correct->Correct, nomatch =',num2str(LEN(5))),'color',COLOR(5,:),'fontsize',8);
text(0,0.6,strcat('Correct->Correct, MATCH =',num2str(LEN(6))),'color',COLOR(6,:),'fontsize',8);
text(0,0.7,strcat('Wrong->Wrong, nomatch =',num2str(LEN(7))),'color',COLOR(7,:),'fontsize',8);
text(0,0.8,strcat('Wrong->Wrong, MATCH =',num2str(LEN(8))),'color',COLOR(8,:),'fontsize',8);



suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12),'-mov'));
cd(Results_dir)
filename = 'PSTH_Simple spike_OCTO_MOV';
print(FF, '-dpdf', filename, '-r400')





%% %%%%%%%%%%% stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_T = Start_T:End_T;
time_M = Start_M:End_M;

% stats8

% WC nomatch
% WC match
% CW nomatch
% CW match
% CC nomatch
% CC match
% WW nomatch
% WW match

CW_DELTA.stats8s = NaN(8,4);
for ii=1:4
    if ii<=2 Psym = Psym_T; time = time_T; end
    if ii>=3 Psym = Psym_M; time = time_M; end
    
    if DELTA.EpochFlag(ii)==1
        ind = find(DELTA.START(ii)<=time & time<=DELTA.END(ii));
        for j=1:8
            CW_DELTA.stats8s(j,ii) = nanmean(Psym(j,ind));
        end
    end
end

   

CW_DELTA.stats8m = NaN(8,4);
for ii=1:4
    if ii<=2 Pmov = Pmov_T; time = time_T; end
    if ii>=3 Pmov = Pmov_M; time = time_M; end
    
    if DELTA.EpochFlag(ii)==1
        ind = find(DELTA.START(ii)<=time & time<=DELTA.END(ii));
        for j=1:8
            CW_DELTA.stats8m(j,ii) = nanmean(Pmov(j,ind));
        end
    end
end












save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');
save(POP_file,'CW_DELTA','-append');





% end
