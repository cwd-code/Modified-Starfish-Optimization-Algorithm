function [Best_Score, BestFit, Convergence_curve] = GWO(SearchAgents,iter_max,lb,ub,dim,fobj)
    % Initialization Boundary
    lb = lb .* ones(1, dim);
    ub = ub .* ones(1, dim);

    % Define Alpha, Beta, Delta Wolf
    Alpha_pos = zeros(1, dim);
    Alpha_score = inf;
    
    Beta_pos = zeros(1, dim);
    Beta_score = inf;
    
    Delta_pos = zeros(1, dim);
    Delta_score = inf;
    
    % Convergence_curve
    Convergence_curve = zeros(1, iter_max);
    
    %% population initialization
    Positions = zeros(SearchAgents, dim);
    fitness = zeros(1, SearchAgents);
    for i = 1:SearchAgents
        Positions(i,:) = (ub - lb) .* rand(1, dim) + lb;
        fitness(i) = fobj(Positions(i,:));
    end
    
    %% Sorting for fitness and finding Alpha, Beta, Delta Wolf
    [SortFitness, indexSort] = sort(fitness);
    Alpha_pos = Positions(indexSort(1),:);
    Alpha_score = SortFitness(1);
    Beta_pos = Positions(indexSort(2),:);
    Beta_score = SortFitness(2);
    Delta_pos = Positions(indexSort(3),:);
    Delta_score = SortFitness(3);
    
    %% Start Iteration
    for t = 1:iter_max
        % Calculate the value of a
        a = 1 + sin(pi/2 + pi * t / iter_max);
        
        for i = 1:SearchAgents
            for j = 1:dim
                %% Updated location based on Alpha Wolf
                r1 = rand();
                r2 = rand();
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_alpha = rand() * sin(rand()) * abs(C1 * Alpha_pos(j) - Positions(i,j));
                X1 = Alpha_pos(j) - A1 * D_alpha;
                
                %% Update location based on Beta Wolf
                r1 = rand();
                r2 = rand();
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_beta = abs(C2 * Beta_pos(j) - Positions(i,j));
                X2 = Beta_pos(j) - A2 * D_beta;
                
                %% Updated location based on Delta Wolf
                r1 = rand();
                r2 = rand();
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_delta = abs(C3 * Delta_pos(j) - Positions(i,j));
                X3 = Delta_pos(j) - A3 * D_delta;
                
                % Update position
                Positions(i,j) = (X1 + X2 + X3) / 3;
            end
            
            % Check boundary
            Positions(i,:) = BoundaryCheck(Positions(i,:), ub, lb, dim);
        end
        
        % Calculation of fitness
        for i = 1:SearchAgents
            fitness(i) = fobj(Positions(i,:));
            
            % Update Alpha, Beta, Delta wolf
            if fitness(i) < Alpha_score
                Alpha_score = fitness(i);
                Alpha_pos = Positions(i,:);
            end
            if fitness(i) > Alpha_score && fitness(i) < Beta_score
                Beta_score = fitness(i);
                Beta_pos = Positions(i,:);
            end
            if fitness(i) > Alpha_score && fitness(i) > Beta_score && fitness(i) < Delta_score
                Delta_score = fitness(i);
                Delta_pos = Positions(i,:);
            end
        end
        
        % Record the current optimal fitness
        Convergence_curve(t) = Alpha_score;
    end
    
    % Return results
    Best_Score = Alpha_score;
    BestFit = Alpha_pos;
end

function X = BoundaryCheck(X, ub, lb, dim)
    for j = 1:dim
        if X(j) > ub(j)
            X(j) = ub(j);
        elseif X(j) < lb(j)
            X(j) = lb(j);
        end
    end
end