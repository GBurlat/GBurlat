% Optimisation Assignment 
% BioComputation
% Gregory Burlat 17006117
% 30/11/2020
% v1.5.0
close all
clear
clc
%%
% Initialising Structures for individuals
parents = struct('gene', {}, 'fitness', 0);
offspring = struct('gene', {}, 'fitness', 0);
tempoff1 = struct('gene', {}, 'fitness', 0);
tempoff2 = struct('gene', {}, 'fitness', 0);

% Initialising variables and arrays
P = 50; % Number of Individuals
N = 10; % Number of values per Gene
mingene = -32.0;
maxgene = 32.0;
generations = 600;
% Mutation Variables
mutationrate = 1/P;
mutstep = 1;

% Arrays for plotting
meanfit = zeros(1,generations);
maxfit = zeros(1,generations);  
minfit = zeros(1,generations);
xaxis = (1:1:generations); % x-axis for plotting

% Loop to randomise parents' genes
for i = 1:P
    for j = 1:N
        parents{i}.gene{j} = mingene + (maxgene - mingene).*rand(1,1);
    end
end

%%
% Generation Loop
w = 1;
while w <= generations
    
    % Calculating fitness of parents and total fitness of parent population
    totalfitp = 0;
    %sigma = 0;
    for i = 1:P
        sigma1 = 0;
        sigma2 = 0;
        for j = 1:N
            x = parents{i}.gene{j};
            sigma1 = sigma1 + x^2;
            sigma2 = sigma2 + cos(2*pi*x);
        end 
        fitness = -20*exp(-0.2*sqrt((1/N)*sigma1))-exp((1/N)*sigma2);
        parents{i}.fitness = fitness;
        totalfitp = totalfitp + fitness;
    end

%%
    % Roulette Selection
    for i = 1:P
        selection_point = totalfitp + (0 - totalfitp).*rand(1,1);
        running_total = 0;
        j = 1;
        while running_total > selection_point
            running_total = running_total + parents{j}.fitness;
            j = j + 1;
        end
        offspring{i} = parents{j-1};            
    end  
    
%%
    % Implementing Crossover
    i = 1;
    while i <= P
        crosspoint = randi([1 N]);
        tempoff1 = offspring{i};
        tempoff2 = offspring{i+1};
        for j = crosspoint:N
            crossbuffer = tempoff1.gene{j};
            tempoff1.gene{j} = tempoff2.gene{j};
            tempoff2.gene{j} = crossbuffer;        
        end
        offspring{i} = tempoff1;
        offspring{i+1} = tempoff2;
        i = i + 2;
    end
%%
    % Mutation
    for i = 1:P
        for j = 1:N
            mutation = rand;
            if mutation <= mutationrate
                alter = 0 + (mutstep-0).*rand(1,1);
                if offspring{i}.gene{j} + alter <= maxgene && offspring{i}.gene{j} - alter >= mingene
                    if randi([0 1]) == 1
                        offspring{i}.gene{j} = offspring{i}.gene{j} + alter;
                    else
                        offspring{i}.gene{j} = offspring{i}.gene{j} - alter;
                    end
                else
                    if offspring{i}.gene{j} + alter >= maxgene
                        offspring{i}.gene{j} = maxgene;  
                    end
                    if offspring{i}.gene{j} - alter <= mingene
                        offspring{i}.gene{j} = mingene;
                    end                        
                end
            end         
        end
    end
 %%
    % Calculating fitness of new mutated offsprings
    totalfitom = 0;
    for i = 1:P
        sigma1 = 0;
        sigma2 = 0;
        for j = 1:N
            x = offspring{i}.gene{j};
            sigma1 = sigma1 + x^2;
            sigma2 = sigma2 + cos(2*pi*x);
        end
         fitness = -20*exp(-0.2*sqrt((1/N)*sigma1))-exp((1/N)*sigma2);
         offspring{i}.fitness = fitness;
         totalfitom = totalfitom + fitness;
    end
    
%%
    % Calculating Mean Fitness
    meanfit(w) = totalfitom/P;

    %Storing Max Fitness of the current Generation
    max = -50;
    maxindex = -50;
    for i = 1:P
        if max < offspring{i}.fitness
           max = offspring{i}.fitness;
           maxindex = i;
        end
    end
    maxfit(w) = max;
    
    %Storing Min Fitness of the current Generation
    min = 0;
    for i = 1:P
        if min > offspring{i}.fitness
           min = offspring{i}.fitness;
        end
    end
    minfit(w) = min;
    
%%
    % Replacing worst fitness candidate of the current population with the
    % best fitness of the previous population
    min = 0;
    minindex = 1;
    for i = 1:P
        if min > parents{i}.fitness
           min = parents{i}.fitness;
           minindex = i;
        end
    end
    
    if offspring{maxindex}.fitness > parents{minindex}.fitness
        offspring{maxindex} = parents{minindex};
    end
    
%%
    % Mutated offspring are the new parents
        parents = offspring;

    % Increase generation loop count
    w = w + 1;
end

% Best fitness individual across all generations
genmin = 0;
for i = 1:w-1
    if genmin > minfit(i)
        genmin = minfit(i);
    end
end

% Best mean fitness population across all generations
minmean = 0;
for i = 1:w-1
    if minmean > meanfit(i)
        minmean = meanfit(i);
    end
end

% Plotting Results
subplot(2,1,1)
plot(xaxis, meanfit, 'g')
xlabel('Generation')
ylabel('Mean Fitness')

subplot(2,1,2)
plot(xaxis, maxfit, 'r')
xlabel('Generation')
ylabel('Max(R) / Min(B) Fitness')

hold on
plot(xaxis, minfit, 'b')