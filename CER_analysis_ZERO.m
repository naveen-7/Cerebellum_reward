

close all;

FLOAG = 4;

% ColorPSTH_n


%% hand movt stuff
if FLOAG==3
    ZERO_manipulandam_NUEVO2
end


%% neural stuff
if FLOAG==1
    ZERO_association_learning_NUEVO
end
if FLOAG==0
    ZERO_Washout_NUEVO
end

if FLOAG==2
    ZERO_motor_learning_NUEVO2
end

if FLOAG==4
    ZERO_association_learning_NUEVO_RT
end



