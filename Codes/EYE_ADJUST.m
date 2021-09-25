
%% adjust X and Y to center the eyegram
% written by naveen at JLG on 12/4/19

function [tempp] = EYE_ADJUST(tempp)



res = 1;


HMEAN = nanmean(tempp);
VMEAN = nanmean(tempp,2);

[val, hpos] = nanmax(HMEAN);
[val, vpos] = nanmax(VMEAN);


figure(); subplot(2,2,1); hold on; plot(HMEAN); plot(VMEAN);


stepy = 50-vpos; stepY = abs(stepy);
if stepy<0
    tempp = tempp(stepY+1*res:end,:);
    tempp = [tempp;zeros(stepY*res,size(tempp,2))];
end

if stepy>0
    tempp = tempp(1:size(tempp,1)-stepY*res,:);
    tempp = [zeros(stepY*res,size(tempp,2)); tempp];
end   


stepx = 50-hpos; stepX = abs(stepx);
if stepx<0
    tempp = tempp(:,stepX*res+1:end);
    tempp = [tempp zeros(size(tempp,1),stepX*res)];
end

if stepx>0
    tempp = tempp(:,1:size(tempp,2)-stepX*res);
    tempp = [zeros(size(tempp,1),stepX*res) tempp];
end  




HMEAN = nanmean(tempp);
VMEAN = nanmean(tempp,2);

subplot(2,2,2); hold on; plot(HMEAN); plot(VMEAN);





delete(gcf);


end