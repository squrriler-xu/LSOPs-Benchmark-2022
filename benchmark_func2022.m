function fit = benchmark_func2022(x, func_num)
% The main function of the 2022 LSOPs benchmark suite
% -------------------------------- Input ----------------------------------
% x: the decision variable
% func_num: the index of the benchmark function.
% -------------------------------- Output ---------------------------------
% fit: the fitness value

global initial_flag

x = x';
if ismember(func_num, [1 : 11, 14, 15])
    fit = Designed_function(x, func_num);
elseif ismember(func_num, [12, 13])
    fit = Nature_function(x, func_num);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BENCHMARK FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = Designed_function(x, func_num)
global initial_flag 
persistent index_base I_all S w b lb ub o R xopt xopt_local

    [D ps] = size(x);
    if (initial_flag == 0)
        filename=sprintf('./cec2022/datafiles_2022/f%02d.mat', func_num);
        load(filename);
        initial_flag = 1;
    end

    idx = checkBounds(x, lb, ub);
    
    bidx = 1;
    comf_set = [];
    for i = 1 : length(S)
        if b(i) == 1
            comf_set = [comf_set, I_all{i}];
        else
            x(I_all{i}, :) = x(I_all{i}, :) - repmat(o{bidx}, 1, ps);
            bidx = bidx + 1;
        end
    end
    comf_set = unique(comf_set);
    x(comf_set, :) = x(comf_set, :) - repmat(o{bidx}, 1, ps);
    
    fit = 0;
    for i = 1 : length(S)
        f = base_functions(x(I_all{i}, :), R{i}, xopt{i}, xopt_local(I_all{i}, :), index_base(i));
        fit = fit + w(i)*f;
    end
    fit(idx) = NaN;
    if ~isempty(idx)
        warning "Warning of boundary crossing.";
    end
end

function fit = Nature_function(x, func_num)
    global initial_flag 
    persistent o lb ub 

    [D ps] = size(x);
    if (initial_flag == 0)
        filename=sprintf('./cec2022/datafiles_2022/f%02d.mat', func_num);
        load(filename);
        initial_flag = 1;
    end

    idx = checkBounds(x, lb, ub);
    x = x-repmat(o, 1, ps);
    
    R = [];
    xopt = [];
    xopt_local = [];
    
    if func_num == 12
        fit = base_functions(x, R, xopt, xopt_local, 13);
    elseif func_num == 13
        fit = base_functions(x, R, xopt, xopt_local, 15);
    end
    fit(idx) = NaN;
    if ~isempty(idx)
        warning "Warning of boundary crossing.";
    end
end


%------------------------------------------------------------------------------
% This function tests a given decision vector against the boundaries of a function.
%------------------------------------------------------------------------------
function indices = checkBounds(x, lb, ub)
    indices = find(sum(x > ub | x < lb) > 0);
end

