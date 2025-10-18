function [fvalbest,xposbest,Curve] = MSFOA(N, T, lb, ub, D, fobj)
    % Starfish optimization algorithm (SFOA) with on-demand distance calculation
    % Modified to compute distances only when needed (reduces time complexity)
    % Parameters:
    %   N: Population size
    %   T: Maximum iterations
    %   D: Problem dimension
    %   fobj: Objective function (assumed O(D) per evaluation)

    %% Initialization
    GP = 0.5;     % Exploration/exploitation parameter 
    k = round(0.05 * N); % Number of worst individuals to reverse
    p_reverse = 0.2; % Probability of reverse learning for others
    if size(ub, 2) == 1
        lb = lb * ones(1, D); 
        ub = ub * ones(1, D);
    end
    fvalbest = inf;
    Curve = zeros(1, T);
    Xpos = rand(N, D) .* (ub - lb) + lb; % Initialize population
    Fitness = zeros(1, N);
    for i = 1:N
        Fitness(i) = feval(fobj, Xpos(i, :));
    end
    [fvalbest, order] = min(Fitness); % Global best
    xposbest = Xpos(order, :);
    newX = zeros(N, D);

    %% Evolution loop
    for t = 1:T
        % Angle calculation (random fluctuation)
        base_amplitude = pi / 2;
        decay_factor = 1 - (t / T);
        theta = base_amplitude * decay_factor * (2 * rand() - 1); 
        theta = max(min(theta, pi), 0);
        tEO = (T - t) / T * cos(theta);

        %% Exploration or Exploitation
        if rand < GP % Exploration
            for i = 1:N
                if D > 5
                    jp1 = randperm(D, 5); % Randomly select 5 dimensions
                    for j = 1:5
                        pm = (2 * rand - 1) * pi;
                        % Compute distance on-demand: Dist(i, order) = norm(Xpos(i,:) - xposbest)
                        dist_i_order = norm(Xpos(i, :) - xposbest);
                        phi = 2.0 * exp(-0.0693 * dist_i_order) + 0.5;
                        origin_X = newX(i, jp1(j));
                        if rand < GP
                            newX(i, jp1(j)) = Xpos(i, jp1(j)) + phi * pm * (xposbest(jp1(j)) - Xpos(i, jp1(j))) * cos(theta);
                        else
                            newX(i, jp1(j)) = Xpos(i, jp1(j)) - phi * pm * (xposbest(jp1(j)) - Xpos(i, jp1(j))) * sin(theta);
                        end
                        w = 0.2 * exp(-5 * t / T) + 0.1; % Dynamic inertia weight
                        % Boundary check
                        if newX(i, jp1(j)) > ub(jp1(j))
                            newX(i, jp1(j)) = w * ub(jp1(j)) + (1 - w) * Xpos(i, jp1(j));
                        elseif newX(i, jp1(j)) < lb(jp1(j))
                            newX(i, jp1(j)) = w * lb(jp1(j)) + (1 - w) * Xpos(i, jp1(j));
                        end
                        newX(i, :) = max(min(newX(i, :), ub), lb); % Full boundary check
                        new_fitness = feval(fobj, newX(i, :));
                        if new_fitness > Fitness(i)
                            newX(i, jp1(j)) = origin_X; % Revert if worse
                        end
                    end
                else % D <= 5
                    jp2 = ceil(D * rand); % Randomly select 1 dimension
                    im = randperm(N, 2); % Select 2 random individuals
                    % Compute distances on-demand: Dist(i, im(1)) and Dist(i, im(2))
                    dist_i_im1 = norm(Xpos(i, :) - Xpos(im(1), :));
                    dist_i_im2 = norm(Xpos(i, :) - Xpos(im(2), :));
                    phi1 = 0.8 * exp(-0.0693 * dist_i_im1) + 0.5;
                    phi2 = 0.8 * exp(-0.0693 * dist_i_im2) + 0.5;
                    rand1 = 2 * rand - 1;
                    rand2 = 2 * rand - 1;
                    origin_X = newX(i, jp2);
                    newX(i, jp2) = tEO * Xpos(i, jp2) + phi1 * rand1 * (Xpos(im(1), jp2) - Xpos(i, jp2)) + phi2 * rand2 * (Xpos(im(2), jp2) - Xpos(i, jp2));
                    w = 0.2 * exp(-5 * t / T) + 0.1;
                    % Boundary check
                    if newX(i, jp2) > ub(jp2)
                        newX(i, jp2) = w * ub(jp2) + (1 - w) * Xpos(i, jp2);
                    elseif newX(i, jp2) < lb(jp2)
                        newX(i, jp2) = w * lb(jp2) + (1 - w) * Xpos(i, jp2);
                    end
                    newX(i, :) = max(min(newX(i, :), ub), lb);
                    new_fitness = feval(fobj, newX(i, :));
                    if new_fitness > Fitness(i)
                        newX(i, jp2) = origin_X;
                    end
                end
            end
        else % Exploitation
            df = randperm(N, 5); % Select 5 random individuals
            dm = zeros(5, D);
            for arm = 1:5
                dm(arm, :) = xposbest - Xpos(df(arm), :); % Direction vectors
            end
            for i = 1:N
                r1 = rand - 0.1; 
                r2 = rand - 0.1;
                kp = randperm(5, 2); % Select 2 random arms
                % Compute distances on-demand: Dist(i, df(kp(1))) and Dist(i, df(kp(2)))
                dist_i_kp1 = norm(Xpos(i, :) - Xpos(df(kp(1)), :));
                dist_i_kp2 = norm(Xpos(i, :) - Xpos(df(kp(2)), :));
                phi1 = 0.8 * exp(-0.0693 * dist_i_kp1) + 0.5;
                phi2 = 0.8 * exp(-0.0693 * dist_i_kp2) + 0.5;
                newX(i, :) = Xpos(i, :) + phi1 * r1 * dm(kp(1), :) + phi2 * r2 * dm(kp(2), :);
                if i == N % Regeneration for last individual
                    newX(i, :) = exp(-t * N / T) .* Xpos(i, :);
                end
                newX(i, :) = max(min(newX(i, :), ub), lb);
            end
        end

        %% Reverse learning mechanism
        [~, idx] = sort(Fitness, 'descend'); % Sort by fitness (worst first)
        worst_idx = idx(1:k); % Worst k individuals
        % Reverse learning for worst individuals (5 attempts per individual)
        for i = worst_idx
            for attempt = 1:5
                dim = randi(D); % Random dimension
                r = rand;
                reverseX = Xpos(i, :);
                reverseX(dim) = ub(dim) + lb(dim) - r * Xpos(i, dim);
                reverse_fitness = feval(fobj, reverseX);
                if reverse_fitness < Fitness(i)
                    newX(i, :) = reverseX;
                end
            end
        end
        % Reverse learning for others (with probability p_reverse)
        for i = 1:N
            if ~ismember(i, worst_idx) && rand < p_reverse
                for attempt = 1:5
                    dim = randi(D);
                    r = rand;
                    reverseX = Xpos(i, :);
                    reverseX(dim) = ub(dim) + lb(dim) - r * Xpos(i, dim);
                    reverse_fitness = feval(fobj, reverseX);
                    if reverse_fitness < Fitness(i)
                        newX(i, :) = reverseX;
                    end
                end
            end
        end

        %% Fitness evaluation and update
        for i = 1:N
            newFit = feval(fobj, newX(i, :));
            if newFit < Fitness(i)
                Fitness(i) = newFit;
                Xpos(i, :) = newX(i, :);
                if newFit < fvalbest
                    fvalbest = newFit;
                    xposbest = Xpos(i, :);
                    order = i; % Update global best index
                end
            end
        end
        Curve(t) = fvalbest;
    end
end