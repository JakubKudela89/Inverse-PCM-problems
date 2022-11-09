function res = Run_DE(Runs,fhd,C,problem_size,funcs,max_nfes,pop_size,optimum)

Rand_Seeds = 10000*rand(1001,1);

Alg_Name=[ 'DE'];

F =0.50*ones(pop_size,problem_size);
Cr=0.90*ones(pop_size,problem_size);

lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];

% fprintf('Running %s algorithm on D= %d\n',Alg_Name, problem_size)

nr_points_to_save = 10;
RecordFEsFactor = round(logspace(log(2*pop_size)/log(10),log(max_nfes)/log(10),nr_points_to_save));


progress = numel(RecordFEsFactor);
val_2_reach = 10^(-8);


for func_no = funcs
%     fprintf('\n-------------------------------------------------------\n')
%     fprintf('Function = %d, Dimension size = %d\n', func_no, problem_size)
    allerrorvals = zeros(progress, Runs);
%     you can use parfor if you have MATLAB Parallel Computing Toolbox
%     parfor run_id = 1 : Runs
    for run_id = 1 : Runs
        run_seed=Rand_Seeds(mod(problem_size*func_no*Runs+run_id-Runs,length(Rand_Seeds)));
        rng(run_seed,'twister');
        Run_RecordFEsFactor=RecordFEsFactor;
        run_funcvals = [];
        run_funccalls = [];
        %% Initialize the main population
        popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
        pop = popold; % the old population becomes the current population
        
        %fitness = feval(fhd,pop',func_no,C);
        fitness = fhd(pop);
        %fitness = fitness';
        nfes = pop_size;
        [bsf_fit_var, minpos] = min(fitness);
        bsf_solution = pop(minpos,:);
        
        run_funcvals = [run_funcvals; bsf_fit_var];
        run_funccalls = [run_funccalls; nfes];
        fprintf('fes: %u, val: %5.3e\n', nfes, bsf_fit_var)
        %% main loop
        while nfes < max_nfes
            pop = popold; % the old population becomes the current population
            R = Gen_R(pop_size,3);
            r1=R(:,2);
            r2=R(:,3);
            r3=R(:,4);
            
            vi = pop(r1, :) + F .* (pop(r2, :) - pop(r3, :));
            
            vi = boundConstraint(vi, pop, lu);
            
            mask = rand(pop_size, problem_size) > Cr; % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);
            
            children_fitness = fhd(ui);
%             children_fitness = feval(fhd, ui', func_no,C);
%             children_fitness = children_fitness';
            
            for i = 1 : pop_size
                nfes = nfes + 1;
                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                    bsf_solution = ui(i,:);
                end
                if nfes > max_nfes; break; end
            end
            
            run_funcvals = [run_funcvals; bsf_fit_var];
            run_funccalls = [run_funccalls; nfes];
            fprintf('fes: %u, val: %5.3e\n', nfes, bsf_fit_var)
            
            
            [fitness, I] = min([fitness, children_fitness], [], 2);
            
            popold = pop;
            popold(I == 2, :) = ui(I == 2, :);
            
        end
        
        if(C(1)==1)
            run_funcvals=run_funcvals-optimum(func_no);
        end
        
        run_funcvals(run_funcvals<val_2_reach)=0;
        
%         fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , run_funcvals(end))
        
        allerrorvals(1:length(run_funcvals), run_id) = run_funcvals;
        allerrorvals(length(run_funcvals):end,run_id) = run_funcvals(end);
        
    end %% end 1 run
    
%     fprintf('min_funvals:\t%e\n',min(allerrorvals(end,:)));
%     fprintf('median_funvals:\t%e\n',median(allerrorvals(end,:)));
%     fprintf('mean_funvals:\t%e\n',mean(allerrorvals(end,:)));
%     fprintf('max_funvals:\t%e\n',max(allerrorvals(end,:)));
%     fprintf('std_funvals:\t%e\n',std(allerrorvals(end,:)));
%     
%     file_name=sprintf('Results\\%s_%s_%s.txt',Alg_Name,int2str(func_no),int2str(problem_size));
    %save(file_name, 'allerrorvals', '-ascii');
    
end %% end 1 function run
res.solution = bsf_solution;
res.bestval = allerrorvals(end);
res.progress_iters = run_funccalls;
res.progress_values = run_funcvals;

end