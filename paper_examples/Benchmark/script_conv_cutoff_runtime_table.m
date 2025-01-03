% Script for comparing the runtime of the solvers (Table 3).

clear all
addpath(genpath(fullfile('../../test_funs/H2004lib/')));
addpath(genpath(fullfile('../../test_funs/TEST29/')));
load('f_min_arr');

fname = 'result-file-name';

disp(['Loading ',fname,'...'])
load(fname);
disp('...done!')

conv_thresh = 10^-4;

%% Create table
T = cell2table(cell(0,6),'VariableNames',{' - ','Gradsamp','DGS','HANSO','SLQPGS','LMBM'});

tmp_cell = cell(2,6);
tmp_cell{1,1} = 'SOGS';
tmp_cell{2,1} = 'Gradsamp';

inds = ([sogs_output.final_f] - f_min_arr' < conv_thresh) & ([gradsamp_output.final_f] - f_min_arr' < conv_thresh);
tmp_cell{1,2} = sum([sogs_output(inds).runtime]);
tmp_cell{2,2} = sum([gradsamp_output(inds).runtime]);

inds = ([sogs_output.final_f] - f_min_arr' < conv_thresh) & ([dgs_output.final_f] - f_min_arr' < conv_thresh);
tmp_cell{1,3} = sum([sogs_output(inds).runtime]);
tmp_cell{2,3} = sum([dgs_output(inds).runtime]);


inds = ([sogs_output.final_f] - f_min_arr' < conv_thresh) & ([hanso_output.final_f] - f_min_arr' < conv_thresh);
tmp_cell{1,4} = sum([sogs_output(inds).runtime]);
tmp_cell{2,4} = sum([hanso_output(inds).runtime]);

inds = ([sogs_output.final_f] - f_min_arr' < conv_thresh) & ([slqpgs_output.final_f] - f_min_arr' < conv_thresh);
tmp_cell{1,5} = sum([sogs_output(inds).runtime]);
tmp_cell{2,5} = sum([slqpgs_output(inds).runtime]);

inds = ([sogs_output.final_f] - f_min_arr' < conv_thresh) & ([lmbm_output.final_f] - f_min_arr' < conv_thresh);
tmp_cell{1,6} = sum([sogs_output(inds).runtime]);
tmp_cell{2,6} = sum([lmbm_output(inds).runtime]);

T = [T;tmp_cell];
disp(T);
    
%% Convert to latex table

fileID = fopen(['table-runtime-',char(datetime("now","Format",'dd-MM-yy_HH-mm-ss')),'-cutoff-',num2str(conv_thresh)],'w');

formatSpec = '%9.1f';
fprintf(fileID,['%% Created via ',pwd,'\\script_conv_cutoff_runtime_table\n']);
fprintf(fileID,['%% Results file: ',fname,'\n']);
fprintf(fileID,['%% conv_thresh = ',num2str(conv_thresh),'\n']);

str = '\\sogs{}';
for i = 2:6
    str = append(str,[' & ',num2str(T{1,i},formatSpec)]);
end
str = append(str,' \\\\ ');
fprintf(fileID,[str,'\n\\hline \n']);

str = 'Other solver';
for i = 2:6
    str = append(str,[' & ',num2str(T{2,i},formatSpec)]);
end
str = append(str,' \\\\ ');
fprintf(fileID,[str,'\n\\hline \n']);