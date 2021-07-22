% Minimisation Function built from "Worksheet_3"
% Biocomputation
% Gregory Burlat 17006117
% 22/11/2020
% v1.0.0
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
P = 500; % Number of Individuals
N = 10; % Number of values per Gene
% f = 10; % Loop constant for fitness calculation
generations = 300;
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
        parents{i}.gene{j} = randi([-512000 512000])/100000;
    end
end

%%
% Generation Loop
w = 1;
while w <= generations
%%
    % Calculating fitness of parents and total fitness of parent population
    totalfitp = 0;
    %sigma = 0;
    for i = 1:P
        sigma = 0;
        for j = 1:N
            x = abs(parents{i}.gene{j});
            sigma = sigma + (x^2)-10*cos(2*pi*x);
        end 
        fitness = 10*N + sigma;
        parents{i}.fitness = fitness;
        totalfitp = totalfitp + fitness;
    end
    
    %%
    % Tournament selection for two parents to generate a better offspring
    totalfito = 0;
    for k = 1:P
        i = randi([1 P]);
        j = randi([1 P]);

        % To avoid choosing the same parent twice
        while i == j
            i = randi([1 P]);
            j = randi([1 P]);
        end

        if parents{i}.fitness > parents{j}.fitness
            offspring{k} = parents{j};
        else
            offspring{k} = parents{i};
        end
        totalfito = totalfito + offspring{k}.fitness;
    end

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

    % Mutation
    for i = 1:P
        for j = 1:N
            mutation = rand;
            if mutation <= mutationrate
                alter = randi([0 mutstep*100])/100;
                if offspring{i}.gene{j} + alter <= 5.12 && offspring{i}.gene{j} - alter >= -5.12
                    if randi([0 1]) == 1
                        offspring{i}.gene{j} = offspring{i}.gene{j} + alter;
                    else
                        offspring{i}.gene{j} = offspring{i}.gene{j} - alter;
                    end
                else
                    if offspring{i}.gene{j} + alter >= 5.12
                        offspring{i}.gene{j} = 5.12;  
                    end
                    if offspring{i}.gene{j} - alter <= -5.12
                        offspring{i}.gene{j} = -5.12;
                    end                        
                end
            end         
        end
    end
    
    % Calculating fitness of new mutated offsprings
    totalfitom = 0;
    for i = 1:P
        sigma = 0;
        for j = 1:N
            x = abs(offspring{i}.gene{j});
            sigma = sigma + (x^2)-10*cos(2*pi*x);
        end
         fitness = 10*N + sigma;
         offspring{i}.fitness = fitness;
         totalfitom = totalfitom + fitness;
    end
        
    % Mutated offspring are the new parents
        parents = offspring;
        
    % Calculating Mean Fitness
    meanfit(w) = totalfitom/P;

    %Storing Max Fitness of the current Generation
    max = 0;
    for i = 1:P
        if max < offspring{i}.fitness
           max = offspring{i}.fitness;
        end
    end
    maxfit(w) = max;
    
    %Storing Min Fitness of the current Generation
    min = 1000;
    for i = 1:P
        if min > offspring{i}.fitness
           min = offspring{i}.fitness;
        end
    end
    minfit(w) = min;

    % Increase generation loop count
    w = w + 1;
end
%%
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