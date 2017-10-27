function [B_Param,B_EYn,B_rXPos,B_rSPos] = yl_COH_filter(COH,Iter)

[row,col] =size(COH);

for i = 1:col
    disp(['Iteration' num2str(i)])
    Yn = COH(:,i);
    n = length(Yn);
    In = ones(n,1);
    valid = ones(n,1);

    %% Set Behavioral Model and Learning Procedure
    %  Create model
    Param = ay_create_state_space(1,0,1,0,1,1,0,1,1);
    Param.Ak = 1;
    Param.Ck = 1;
    % Set learning parameters

    Param = ay_set_learning_param(Param,Iter,0,1,1,1,1,0,1,2,1);
    % EM
    [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb]=ay_em([2 0],[],In,0,Yn,[],Param,valid);

    % Extract estimating data
    B_XSmt{i} = rXSmt;
    B_SSmt{i} = rSSmt;

    B_rXPos{i} = rXPos;
    B_rSPos{i} = rSPos;

    B_Param{i} = Param;
    B_EYn{i} = EYn;
    B_rYn{i} = rYn;
end



end