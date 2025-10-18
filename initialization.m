
%% Population initialization function
function X = initialization(pop,dim,ub,lb)
    %pop:population size
    %dim:Each individual dimension
    %ub:Upper bound of the variable for each dimension with dimensions [1,dim];
    %lb:Lower bound of the variable for each dimension with dimensions [1,dim]; 
    %X:Exported populations£¬with dimensions [pop,dim];
    X = zeros(pop,dim); %Pre-allocate space for X
    for i = 1:pop
       for j = 1:dim
           X(i,j) = (ub(j) - lb(j))*rand() + lb(j);  %Generate a random number between [lb,ub]
       end
    end
end
