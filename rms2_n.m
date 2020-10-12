function rms = rms2_n(X,Y)

if size(X)==size(Y)
    for i=1:length(X)
        RMS(i) = (X(i)-Y(i))^2;
    end
end

rms = sqrt(nansum(RMS));


end