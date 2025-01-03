% A function for displaying the performance of all solvers as a table

function perf_table(best_fval_overall,sogs_output,gradsamp_output,dgs_output,hanso_output,slqpgs_output,lmbm_output)

T = table({'SOGS';'Gradsamp';'DGS';'Hanso';'SLQPGS';'LMBM'},...
    [ ...
    sogs_output.f_eval; ...
    gradsamp_output.f_eval; ...
    dgs_output.f_eval; ...
    hanso_output.f_eval; ...
    slqpgs_output.f_eval; ...
    lmbm_output.f_eval ...
    ],...
    [ ...
    sogs_output.subgrad_eval; ...
    gradsamp_output.subgrad_eval; ...
    dgs_output.subgrad_eval; ...
    hanso_output.subgrad_eval; ...
    slqpgs_output.subgrad_eval; ...
    lmbm_output.subgrad_eval ...
    ],...
    [ ...
    sogs_output.subhess_eval; ...
    gradsamp_output.subhess_eval; ...
    dgs_output.subhess_eval; ...
    hanso_output.subhess_eval; ...
    slqpgs_output.subhess_eval; ...
    lmbm_output.subhess_eval ...
    ],...
    round([ ...
    sogs_output.final_f; ...
    gradsamp_output.final_f; ...
    dgs_output.final_f; ...
    hanso_output.final_f; ...
    slqpgs_output.final_f; ...
    lmbm_output.final_f ...
    ],4,"significant"),...
    [ ...
    sogs_output.final_f - best_fval_overall; ...
    gradsamp_output.final_f - best_fval_overall; ...
    dgs_output.final_f - best_fval_overall; ...
    hanso_output.final_f - best_fval_overall; ...
    slqpgs_output.final_f - best_fval_overall; ...
    lmbm_output.final_f - best_fval_overall ...
    ],...
    [ ...
    sogs_output.runtime; ...
    gradsamp_output.runtime; ...
    dgs_output.runtime; ...
    hanso_output.runtime; ...
    slqpgs_output.runtime; ...
    lmbm_output.runtime ...
    ]);
T.Properties.VariableNames = {'Algorithm','f','subgrad_f','subhess_f','Final value','Difference','Time'};
disp(T);

end
