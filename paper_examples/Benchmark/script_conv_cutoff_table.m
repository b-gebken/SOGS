% Script for comparing the number of evaluations of all solvers (Table 1
% and Table 2). 

clear all
addpath(genpath(fullfile('../../test_funs/H2004lib/')));
addpath(genpath(fullfile('../../test_funs/TEST29/')));
load('f_min_arr');

fname = 'result-file-name';

disp(['Loading ',fname,'...'])
load(fname);
disp('...done!')

conv_thresh = 0; % For Table 1
% conv_thresh = 10^-4; % For Table 2

%% Create table
T = cell2table(cell(0,14),'VariableNames',{'Prob. No.','Winner',...
    'eval SOGS','acc sogs',...
    'eval Gradsamp','acc GS',...
    'eval DGS','acc DGS',...
    'eval Hanso','acc Hanso',...
    'eval SLQPGS','acc SLQPGS',...
    'eval LMBM','acc LMBM'});

for i = 1:numel(problem_arr)
    disp(['i = ',num2str(i),'...']);

    sogs_fx_arr = compute_fx_arr(problem_arr(i).f,sogs_output(i).x_arr);
    sogs_ind = min(find(sogs_fx_arr - f_min_arr(i) >= conv_thresh,1,'last')+1,numel(sogs_fx_arr));
    if(isempty(sogs_ind))
        sogs_ind = 1;
    end
    sogs_subgrad_arr = tril(ones(sogs_output(i).num_iter+1)) * [sogs_output(i).numsample_cell{:}]'; 

    gradsamp_f_eval_pts = compute_fx_arr(problem_arr(i).f,gradsamp_output(i).eval_pts);
    dgs_f_eval_pts = compute_fx_arr(problem_arr(i).f,dgs_output(i).eval_pts);
    hanso_f_eval_pts = compute_fx_arr(problem_arr(i).f,hanso_output(i).eval_pts);
    slqpgs_f_eval_pts = compute_fx_arr(problem_arr(i).f,slqpgs_output(i).eval_pts);
    lmbm_f_eval_pts = compute_fx_arr(problem_arr(i).f,lmbm_output(i).eval_pts);
    
    gradsamp_ind = find(gradsamp_f_eval_pts - f_min_arr(i) <= conv_thresh,1,'first');
    dgs_ind = find(dgs_f_eval_pts - f_min_arr(i) <= conv_thresh,1,'first');
    hanso_ind = find(hanso_f_eval_pts - f_min_arr(i) <= conv_thresh,1,'first');
    slqpgs_ind = find(slqpgs_f_eval_pts - f_min_arr(i) <= conv_thresh,1,'first');
    lmbm_ind = find(lmbm_f_eval_pts - f_min_arr(i) <= conv_thresh,1,'first');
    
    if(isempty(gradsamp_ind))
        gradsamp_ind = numel(gradsamp_f_eval_pts);
    end
    if(isempty(dgs_ind))
        dgs_ind = numel(dgs_f_eval_pts);
    end
    if(isempty(hanso_ind))
        hanso_ind = numel(hanso_f_eval_pts);
    end
    if(isempty(slqpgs_ind))
        slqpgs_ind = numel(slqpgs_f_eval_pts);
    end
    if(isempty(lmbm_ind))
        lmbm_ind = numel(lmbm_f_eval_pts);
    end

    tmp_cell = {sogs_subgrad_arr(sogs_ind);...
        gradsamp_ind;...
        dgs_ind;...
        hanso_ind;...
        slqpgs_ind; ...
        lmbm_ind};
    
    acc_arr = [sogs_fx_arr(sogs_ind) - f_min_arr(i); ...
        gradsamp_f_eval_pts(gradsamp_ind) - f_min_arr(i); ...
        dgs_f_eval_pts(dgs_ind) - f_min_arr(i); ...
        hanso_f_eval_pts(hanso_ind) - f_min_arr(i); ...
        slqpgs_f_eval_pts(slqpgs_ind) - f_min_arr(i); ...
        lmbm_f_eval_pts(lmbm_ind) - f_min_arr(i)];
    
    eval_arr = cell2mat(tmp_cell);
    eval_arr(acc_arr >= conv_thresh) = Inf;
    
    min_numb = min(eval_arr,[],1);
    if(isinf(min_numb))
        winner_ind = {[]};
    else
        winner_ind = {find(eval_arr == min_numb)'};
    end
    
    new_row_cell = {i,winner_ind,sogs_subgrad_arr(sogs_ind),sogs_fx_arr(sogs_ind) - f_min_arr(i),...
        gradsamp_ind,gradsamp_f_eval_pts(gradsamp_ind) - f_min_arr(i),...
        dgs_ind,dgs_f_eval_pts(dgs_ind) - f_min_arr(i),...
        hanso_ind,hanso_f_eval_pts(hanso_ind) - f_min_arr(i),...
        slqpgs_ind,slqpgs_f_eval_pts(slqpgs_ind) - f_min_arr(i),...
        lmbm_ind,lmbm_f_eval_pts(lmbm_ind) - f_min_arr(i)};
    
    T = [T;new_row_cell];
    
end

disp(T);
    
%% Convert to latex table
T_mat = table2array(T(:,[1,3:end]));

fileID = fopen(['table-',char(datetime("now","Format",'dd-MM-yy_HH-mm-ss')),'-cutoff-',num2str(conv_thresh)],'w');

formatSpec = '%9.1e';
fprintf(fileID,['%% Created via ',pwd,'\\script_conv_cutoff_table\n']);
fprintf(fileID,['%% Results file: ',fname,'\n']);
fprintf(fileID,['%% conv_thresh = ',num2str(conv_thresh),'\n']);
for i = 1:numel(problem_arr)
    str = [num2str(i),'. & '];
    for j = 1:6
        if(any(j == T.Winner{i}) && conv_thresh ~= 0)
            str = append(str,['\\textbf{',num2str(T_mat(i,2 + 2*(j-1))),'} & \\textbf{',num2str(T_mat(i,2 + 2*j - 1),formatSpec),'}']);
        elseif(T_mat(i,2 + 2*j - 1) >= conv_thresh && conv_thresh ~= 0)
            str = append(str,'- & -');
        else
            str = append(str,[num2str(T_mat(i,2 + 2*(j-1))),' & ',num2str(T_mat(i,2 + 2*j - 1),formatSpec)]);
        end
        
        if(j == 6)
            str = append(str,' \\\\ ');
        else
            str = append(str,' & ');
        end
    end
    
    str = replace(str,'0e+00','0');
    fprintf(fileID,[str,'\n\\hline \n']);
end

%% Functions

function fx_arr = compute_fx_arr(f,x_arr)
    fx_arr = zeros(1,size(x_arr,2));
    for j = 1:size(x_arr,2)
        fx_arr(j) = f(x_arr(:,j));
    end
end