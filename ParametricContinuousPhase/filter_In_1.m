%% Load data
load('subset_COH.mat')

[row,col] =size(COH)

%% Extract data at each frequency
% Store 200 samples in each cell of Y_COH
Y_COH = cell(col,1);
for i = 1:col
    Y_COH{i} = COH(:,i)';
end


for i = 1:col
disp(['Iteration' num2str(i)])
Yn = Y_COH{i};
n = length(Yn);
In = ones(n,1);
valid = ones(n,1);

%% Set Behavioral Model and Learning Procedure
%  Create model
Param = ay_create_state_space(1,0,1,0,1,1,0,1,1);
Param.Ak = 1;
% Set learning parameters
Iter = 30;
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