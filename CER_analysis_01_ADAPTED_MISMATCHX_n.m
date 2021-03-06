%% Takes data from REX and Spike2 and merges them and makes the data ready
%% for further analysis

% Created by Anna and Adam
% Modified by Naveen on 11/30/15 at CUMC
% Adapted by Naveen on 12/29/16 at CUMC
% Xd by Naveen on 08/11/17 at CUMC
% Edited on 8/11/19

% function CER_analysis_01_ADAPTED_MISMATCHX_n


% clc;
% clear all;
% close all;


%% Setup directories--------------------------------------------------------

disp('!!!!!  CER_analysis_01_ADAPTED_MISMATCHX_n has started running  !!!!!')


% CODES ----------------------------------------------------

START_ = 800;
FIXON_ = 1003;
TARG_ = 1007;
CORR_ = 1012;
WRONG_ = 1013;
EOT_ = 1014;
STARTSPKREC_ = 1095;
ENDSPKREC_ = 1096;
END_ = 801;




%%
%% ONE: ALIGN REX AND SPIKE2
disp('*************************************************')
disp('*************************************************')
disp('ALIGNING REX AND SPIKE2')

data={};        % final data structure
% column 1: unsorted spike times
% column 2: event times
% column 3: sorted spike times
% column 4: rex sorted
Rexdiff=[];     % time differences between 1095 and 1096 for each trial from REX
Spikediff=[];   % time differences between 1095 and 1096 for each trial from Spike2
diffdiff=[];    % the difference between the two previously mentioned differences



% DO NOT TAKE FIRST TRIAL
% find time of first 1001 eCode to calculate start time of first trial
startTimes=[NaN];
endTimes=[NaN];
loc95=NaN;
for i=size(Trials(1).Events,2):-1:1
    if Trials(1).Events(i).Code == START_
        loc95=i+1;
        startTimes=[double(Trials(1).Events(i).Time)];
        endTimes=[double(Trials(1).Events(size(Trials(1).Events,2)).Time)];
        Rexdiff=[double(Trials(1).Events(size(Trials(1).Events,2)-1).Time)-double(Trials(1).Events(i+1).Time)];
        break;
    end
end

for i=2:size(Trials,2)
    startTimes=[startTimes; double(Trials(i).Events(1).Time)];
    endTimes=[endTimes; double(Trials(i).Events(size(Trials(i).Events,2)).Time)];
    Rexdiff=[Rexdiff; double(Trials(i).Events(size(Trials(i).Events,2)-1).Time)-double(Trials(i).Events(2).Time)];
end

for i=2:2:size(Trigg.times,1)
    Spikediff=[Spikediff; Trigg.times(i)-Trigg.times(i-1)];
end
Spikediff=Spikediff*1000;



% check if we have collected the wrong intervals from Spike
count = 0;
start3=0;
for i=1:size(Spikediff,1)
    if Spikediff(i) > 779 && Spikediff(i) < 801
        count = count+1;
    else
        count = 0;
    end
    if count >= 5
        Spikediff=[];
        for i=3:2:size(Trigg.times,1)
            Spikediff=[Spikediff; Trigg.times(i)-Trigg.times(i-1)];
        end
        Spikediff=Spikediff*1000;
        start3=1;
        break;
    end
end




sBig=[];
if size(Spikediff,1) > size(Rexdiff,1)
    sBig=1;
    A=Spikediff;
    B=Rexdiff;
else
    sBig=0;
    B=Spikediff;
    A=Rexdiff;
end

%find trial offset between REX and Spike
exit = 0;
match=0;
offset=0; % how far the shorter array of A/B is shifted to match the longer array
while match==0 && exit == 0
    diffdiff=(A(offset+1:offset+size(B,1))-B);
    for i=1:20
        if diffdiff(i)>1
            offset=offset+1;
            if offset+size(B,1)>size(A,1)
                exit = 1;
            end
            break;
        end
        if i==20
            match=1;
        end
    end
end
% % 
% % if match==0
% %     offset=-1;
% %     while match==0
% %         diffdiff=A(1:size(B,1)+offset)-B(1-offset:size(B,1));
% %         for i=1:20
% %             if diffdiff(i)>1
% %                 offset=offset-1;
% %                 break;
% %             end
% %             if i==20
% %                 match=1;
% %             end
% %         end
% %     end
% % end


timeOffset=[]; % Rex-1000*Spike.  Time from 1000*Spike to Rex
if sBig && offset < 0
    timeOffset=double(Trials(1-offset).Events(2).Time)-Trigg.times(1+start3)*1000;
elseif sBig && offset >= 0
    timeOffset=double(Trials(1).Events(loc95).Time)-Trigg.times(2*offset+1+start3)*1000;
elseif ~sBig && offset >= 0
    timeOffset=double(Trials(1+offset).Events(2).Time)-Trigg.times(1+start3)*1000;
elseif ~sBig && offset < 0
    timeOffset=double(Trials(1).Events(loc95).Time)-Trigg.times(2*(-offset)+1+start3)*1000;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TWO: SPIKE SORTING
disp('*************************************************')
disp('*************************************************')
disp('SIMPLE SPIKE SORTING')

%
% UNIQUE_CODES = unique(CHAN.codes(:,1));
% max_sspikes   = max(UNIQUE_CODES);


% Repetitive Iteration or Solo Process Algorithm--------------------



if  isempty(CHAN.codes)
    disp(strcat('------There is no Spike Code available--------'));
    disp('>>> Proceeding >>>');
    CODE_S = 0;
    max_sspikes   = 0;
    sUNIQUE_CODES = 0;
    types=0;
else
    UNIQUE_CODES = unique(CHAN.codes(:,1));
    types=max(CHAN.codes(:,1));
    
    if max(UNIQUE_CODES)==1
        disp(strcat('------There is only ONE Spike Code available--------'));
        disp('>>> Proceeding >>>');
        CODE_S = 1;
        max_sspikes   = max(UNIQUE_CODES);
        
    else
        FlagF=0;
        fprintf('----------There are  %d  codes available---------\n',max(UNIQUE_CODES));
        fprintf(strcat('-----------  in this ', String,' spike file  ---------\n'));
        
        max_sspikes   = max(UNIQUE_CODES);
        sUNIQUE_CODES = UNIQUE_CODES;
        
        while FlagF==0
            disp('----------Please enter the code to analyze---------');
            CODE_S = str2num(input('Here: ','s'));
            
            if CODE_S==0
                break;
            else
                
                for kkk=1:length(UNIQUE_CODES)
                    
                    if CODE_S == UNIQUE_CODES(kkk)
                        FillArray(kkk)=1;
                    else
                        FillArray(kkk)=0;
                    end
                    
                end
                
                
                if sum(FillArray)>0
                    FlagF=1;
                    break;
                    
                else
                    FlagF=0;
                end
                
            end
        end
        
    end
    
    
    
end





trial=1;
spikes=[];
all_shape=[];

sorted=cell(types,1);
initial_shape=cell(types,1);
SHAPES={};
i=1;
while i<size(CHAN.times,1)
    if CHAN.times(i)*1000+timeOffset<startTimes(1)
        i=i+1;
    else
        break;
    end
end

while i<size(CHAN.times,1)
    %     trial
    if (trial+1)>size(startTimes,1) || CHAN.times(i)*1000+timeOffset<startTimes(trial+1)
        spikes=[spikes;CHAN.times(i)*1000+timeOffset-startTimes(trial)];
        all_shape=[all_shape;CHAN.values(i,:)];
        if CHAN.codes(i,1) ~= 0
            sorted{CHAN.codes(i,1)}=[sorted{CHAN.codes(i,1)};CHAN.times(i)*1000+timeOffset-startTimes(trial)];
            initial_shape{CHAN.codes(i,1),:}=[initial_shape{CHAN.codes(i,1)};CHAN.values(i,:)];
        end
    else
        data{trial,1}=spikes;
        data{trial,3}=sorted;
        SHAPES{trial,1}=all_shape;
        SHAPES{trial,2}=initial_shape;
        trial=trial+1;
        spikes=[];
        all_shape=[];
        sorted=cell(types,1);
        initial_shape=cell(types,1);
        continue;
    end
    i=i+1;
end
data{trial,1}=spikes;
data{trial,3}=sorted;
SHAPES{trial,1}=all_shape;
SHAPES{trial,2}=initial_shape;


% add Event times
for i=1:size(Trials,2)
    data{i,2}=Trials(i).Events;
    for j=1:size(Trials(i).Events,2);
        data{i,2}(j).Time=data{i,2}(j).Time-startTimes(i);
    end
    if isempty(data{i,3})
        data{i,3}={[]};
        SHAPES{i,2}={[]};
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  THREE Getting all the events
disp('*************************************************')
disp('*************************************************')
disp('GETTING ALL THE EVENTS')

for k=1:size(data,1)
    datcode=data{k,2};
    sp=data{k,3};
    %     datsp=round(sp{1})
    %     [the_struct(xx).Events.Code; the_struct(xx).Events.Time]
    dat(k).Events.Code=[datcode.Code];
    dat(k).Events.Time=double([datcode.Time]);
    
    for y=1:length(data{k,3})
        
        dat(k).Units(y).Code=y;
        datsp=round(sp{y});
        dat(k).Units(y).Times=double(round(datsp))';
    end
    eventpack= [double([datcode.Code])' double([datcode.Time])'];
end








% filenm=[filename '_SSspk']
leng = length(Trials); % the # trials without trial 1.
for xx = 1:leng
    
    
     xx
    
    event_packet = [double([Trials(xx).Events.Code]); double([Trials(xx).Events.Time])]';
    event_packet(find(event_packet<0))=NaN;
    
    num_events = size(event_packet,1);
    if  num_events <=50
        for yy = 2:num_events
            Trials(xx).Events(yy).Time=(Trials(xx).Events(yy).Time-Trials(xx).Events(1).Time);
        end
        
        event_packets=[[Trials(xx).Events.Code]' [Trials(xx).Events.Time]' ];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_units = length(Trials(xx).Units);%.Times);
        %     if Trials(xx).Units.Times(1)>10000
        for yy = 1:num_units
            Trials(xx).Units(yy).Times=(Trials(xx).Units(yy).Times-Trials(xx).Events(1).Time);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Put in array number
        if ~isempty(event_packet(find((event_packet(:,1) >= 2000) & (event_packet(:,1) < 3000)),1) - 2000)
            Trials(xx).Array = event_packet(find((event_packet(:,1) >= 2000) & (event_packet(:,1) < 3000)),1) - 2000;
        else
            Trials(xx).Array=NaN;
        end
        
        if ~isempty (Trials(xx).Events(find(event_packet(:,1) == 1020)))
            Trials(xx).barrelease=1;
        end
        % quewto loop "if" ho dovuto aggiungerlo perche'
        % molti trials non c'era il code 1007?????????
        % 16 aprile 02
        if ~isempty (Trials(xx).Events(find(event_packet(:,1) == 1007)))
            rf_on_time = double(event_packet(find(event_packet(:,1)==1007),2));
            if size(rf_on_time,1)>1
                rf_on_time=rf_on_time(1);
            end
            %             if rf_on_time > 0
            
            % Fix y-eye position by shift_amount
            % This ASSUMES that the eye postion (eyev) is signal # 2.
            if ~isempty(Trials(xx).Signals)
                temp = [double(Trials(xx).Signals(2).Signal)]';
                temp = temp;% + shift_amount;
                for yy = 1:(length(temp))
                    Trials(xx).Signals(2).Signal(yy) = temp(yy);
                end
            end
            
            % Now find RF object information for each trial
            first_patcd = min(find((event_packet(:,1)>3000)& (event_packet(:,1)<4000) &(event_packet(:,2)==rf_on_time))); % The 3000 code gives the # of objects following
            num_objs = double(event_packet(first_patcd,1)) - 3000;
            num_objs=1;
            last_patcd = max(find((event_packet(:,2) == rf_on_time)&(event_packet(:,1) > 5000))); % this is to double check
            
            if (mod(((last_patcd + 1) - (first_patcd + 1)),6)>0)
                %err_text = ['In Trial ' num2str(xx) 'There are ' num2str((last_patcd + 1) - first_patcd) ' objects.'];
                %errordlg('There are not an even number of codes','Warning');
            else
                if ((((last_patcd + 1) - (first_patcd + 1))/6)~= num_objs)
                    %  errordlg('There are a different number of rf_objects to what is specificed','Warning');
                end
            end
            


            counter = find((event_packet(:,1)>7000)& (event_packet(:,1)<8000));

            if event_packet(counter,1) - 7000 >0
                Trials(xx).pattern=event_packet(counter,1) - 7000;
            else
                Trials(xx).pattern=NaN;
            end
            
            Trials(xx).rf_on_time = double(event_packet(find(event_packet(:,1)==1007),2));
            if size(Trials(xx).rf_on_time,1)>1
                Trials(xx).rf_on_time= Trials(xx).rf_on_time(1);
            end
            if ~isempty(event_packet(find(event_packet(:,1)==1010),2))
                Trials(xx).rf_off_time = double(event_packet(find(event_packet(:,1)==1010),2));
            else
                Trials(xx).rf_off_time = NaN;
            end
            
            
            if ~isempty(event_packet(find(event_packet(:,1)==1093),2))
                Trials(xx).abort=double(event_packet(find(event_packet(:,1)==1093),2));
            else
                Trials(xx).abort=NaN;
            end
            
            
            if ~isempty(event_packet(find(event_packet(:,1)==1003),2))
                Trials(xx).FP_on=double(event_packet(find(event_packet(:,1)==1003),2));
            else
                Trials(xx).FP_on=NaN;
            end
            
            if ~isempty(event_packet(find(event_packet(:,1)==1012),2))
                Trials(xx).correct = 1;
                Trials(xx).reward= double(event_packet(find(event_packet(:,1)==1012),2));%-double([Trials(xx).Events(1).Time]);
                DUM = Trials(xx).reward-Trials(xx).rf_on_time;
                DUM = DUM(find(DUM>0));
                if isempty(DUM) DUM = NaN; else 
                    if length(DUM)>1 DUM = min(DUM); end 
                end
                Trials(xx).MRT= DUM; %%%%%%%%%%%
            else
                Trials(xx).correct = NaN;
                Trials(xx).reward= NaN;
                Trials(xx).MRT=NaN;
            end
            
            if ~isempty(event_packet(find(event_packet(:,1)==1013),2))
                Trials(xx).correct = 0;
                Trials(xx).reward=NaN;
                error_time=double(event_packet(find(event_packet(:,1)==1013),2));
                DUM = error_time-Trials(xx).rf_on_time;
                DUM = DUM(find(DUM>0)); 
                
                if isempty(DUM) DUM = NaN; else 
                    if length(DUM)>1 DUM = min(DUM); end 
                end
                Trials(xx).MRT= DUM; %%%%%%%%%%%
            end
            
        else
            Trials(xx).pattern=NaN;
            Trials(xx).rf_on_time =NaN;
            Trials(xx).rf_off_time = NaN;
            Trials(xx).abort=NaN;
            Trials(xx).FP_on=NaN;
            Trials(xx).correct = NaN;
            Trials(xx).reward= NaN;
            Trials(xx).MRT=NaN;
            Trials(xx).obj.pattern = NaN;
            Trials(xx).obj.xpos = NaN;
            Trials(xx).obj.ypos = NaN;
            Trials(xx).obj.red = NaN;
            Trials(xx).obj.green = NaN;
            Trials(xx).obj.blue = NaN;
            Trials(xx).Array=NaN;
            
        end
    else
        Trials(xx).pattern=NaN;
        Trials(xx).rf_on_time =NaN;
        Trials(xx).rf_off_time = NaN;
        Trials(xx).abort=NaN;
        Trials(xx).FP_on=NaN;
        Trials(xx).correct = NaN;
        Trials(xx).reward= NaN;
        Trials(xx).MRT=NaN;
        Trials(xx).obj.pattern = NaN;
        Trials(xx).obj.xpos = NaN;
        Trials(xx).obj.ypos = NaN;
        Trials(xx).obj.red = NaN;
        Trials(xx).obj.green = NaN;
        Trials(xx).obj.blue = NaN;
        Trials(xx).Array=NaN;
        
    end
    
   
    
end













ar=[Trials.Array]';
mrt=[Trials.MRT]';
ptr=[Trials.pattern]';

thistr=(~isnan(ar(:,1)) & ~isnan(mrt(:,1)) & ~isnan(ptr(:,1)));
Trials=Trials(thistr);
data=data(thistr,:);
SHAPES=SHAPES(thistr,:);
dat=dat(thistr);


%% in some trials there is not the Units Code; take out those trials
wr=1;
for u=1:length(Trials);
    if ~isempty(Trials(u).Units);
        Trials(wr)=Trials(u);
        wr=wr+1;
    end
    allTrials(u,1)=Trials(u);
end




for k = 1:length(allTrials)
    temp = allTrials(k).Events;
    val=0;
    for j=1:length(temp)
        if temp(j).Code==1012
            val=1;
        end
        if temp(j).Code==1013
            val=2;
        end
    end
    
    OUTCOME(k,1) = val;
    ARRAY(k,1) = Trials(k).Array;
    PATTERN(k,1) = Trials(k).pattern;
    
end


data_CH = data;
shapes_CH = SHAPES;



%% SIX Some conversions -------
%  fn=[filenm '_spk']
% len=length(FileName1)-4;
% filenm = strcat('Data_',FileName1(1:len));
% cd(PathName1);

clear Infos
Infos = NaN(size(data_CH,1),20);
for i=1:size(dat,2)
    %     i
    clear Temp;
    Temp = [[double(dat(i).Events.Code)];[double(dat(i).Events.Time)]]';
    
    if ~isempty(Temp(find((Temp(:,1)==800))))   & ~isempty(Temp(find((Temp(:,1)==1003)))) & ...
            ~isempty(Temp(find((Temp(:,1)==1007))))  & ~isempty(Temp(find((Temp(:,1)==1014)))) & ...
            ~isempty(Temp(find((Temp(:,1)==801))))   & ...
            (~isempty(Temp(find((Temp(:,1)==1012)))) | ~isempty(Temp(find((Temp(:,1)==1013)))))
        %
        
        
        Infos(i,1)  = Temp(find(Temp(:,1)== 800),2); %Start
%         Infos(i,2)  = Temp(find(Temp(:,1)==1095),2); %Start Spike Recording
        tempp = Temp(find(Temp(:,1)==1003),2);
        Infos(i,3)  = tempp(tempp>0);  %Monkey starts to fixate
        if size(Temp(find(Temp(:,1)==1007),2),1)>1
            bugvar=Temp(find(Temp(:,1)==1007),2);
            Infos(i,4)  = bugvar(1);
        else
            Infos(i,4)  = Temp(find(Temp(:,1)==1007),2); %Target
        end
        Infos(i,5)  = PATTERN(i); %Target shape
        Infos(i,9)  = ARRAY(i); %Hand (0-Right 1-Left)
        Infos(i,10) = OUTCOME(i);    %Correct = 1; Wrong = 2
        tempp = Temp(find((Temp(:,1)==1012) | (Temp(:,1)==1013)),2);
        tempp = tempp(tempp>0 & tempp<6000);
        if isempty(tempp) tempp=NaN; end
        Infos(i,11) = tempp;  %Time of correct or wrong
        % % %         Infos(i,11) = Temp(find(Temp(:,1)==1014),2)-1;  %Time of correct or wrong
        %         Infos(i,12) = Temp(find(Temp(:,1)==1014),2); %End of trial
        tempp = Temp(find(Temp(:,1)== 801),2);
%         Infos(i,13) = tempp(tempp>0); %End of recording
        Infos(i,14) = Infos(i,11)-Infos(i,4);    %RT
        Infos(i,15) = Infos(i,13)-Infos(i,1);    %Trial duration
        Infos(i,16) = 0;    %Category
        
        % % %          clear TEMP;
        % % %         if ~isempty(Temp(find(Temp(:,1)==1050),2)) & ~isempty(Temp(find(Temp(:,1)==1051),2))
        % % %             TEMP = Temp(find(Temp(:,1)==1051),2)-Temp(find(Temp(:,1)==1050),2);
        % % %             Infos(i,17) =  TEMP(find(TEMP>0));
        % % %
        % % %             Infos(i,11) = Temp(find(Temp(:,1)==1050),2);  %Time of correct or wrong
        % % %         end
        
    else
        
        
        data_CH{i,1} = [];
        data_CH{i,3}{1,1} = [];
        shapes_CH{i,1} = [];
        shapes_CH{i,2}{1,1} = [];
        
        
    end
end





disp('!!! END OF CODE !!!');
disp('!!! Data ready for further analysis !!!');


CODE1_Pass = 1;




% end
