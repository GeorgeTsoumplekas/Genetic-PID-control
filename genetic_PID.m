function [chromosome, min_fit] = genetic_PID()
%function implementing the genetic algorithm

    %number of individuals in the population
    N = 100;
    
    %number of selected individuals
    M = 20;
    
    %number of elite children
    elite_count = 10;
    
    %maximum number of generations (set to 100*#variables)
    max_gens = 600;
    
    %maximum number of generations with no change in their best individual
    max_stall_gens = 20;
    
    %max time of observations (in s)
    tf = 8;
    
    %matrices containing the chromosomes of the parent generation, of the
    %selected parents and of the children generation
    parent_popV = zeros(N,6);
    children_popV = zeros(N,6);
    select_popV = zeros(M,6);
    
    %1st step: create the initial population
    parent_popV(:,1:2) = unifrnd(0,100,N,2);
    parent_popV(:,3:4) = unifrnd(0,10,N,2);
    parent_popV(:,5:6) = unifrnd(0,100,N,2);
    
    fit_val = zeros(N,1);
    
    for i = 1:N
        [xV,tV,uV] = simclosedloop(parent_popV(i,1),parent_popV(i,2),parent_popV(i,3),parent_popV(i,4),parent_popV(i,5),parent_popV(i,6),tf);
        
        x1V = xV(:,1);
        x3V = xV(:,3);
        u1V = uV(:,1);
        u2V = uV(:,2);
        
        fit_val(i) = fitness_function(x1V,x3V,u1V,u2V,tV);
    end
    
    curr_gen = 0;
    curr_stall_gen = 0;
    
    %chromosomes with the best fitness value
    elit_par = min(fit_val);
    elit_child = 0;
    
    %main loop of the algorithm
    while curr_gen < max_gens && curr_stall_gen < max_stall_gens
    
        curr_gen = curr_gen + 1;
        
        if curr_gen ~= 1
            parent_popV = children_popV;
            elit_par = elit_child;
        end
        
        %elite individuals pass to the next generation right away
        [~,orig_indx] = sort(fit_val);
        elite_indx = orig_indx(1:elite_count);
        
        children_popV(1:elite_count,:) = parent_popV(elite_indx,:);
    
    
        %selection process applying the simulation of a roulette to choose
        %the points
        fit_min = min(fit_val);
        fit_max = max(fit_val);
    
        scaled_fit = 1 - (fit_val-fit_min)/(fit_max-fit_min);

        fit_tot = sum(scaled_fit);

        select_area = zeros(N,1);

        select_area(1) = 0;
        for i = 2:N
           select_area(i) =  select_area(i-1) + scaled_fit(i-1)/fit_tot;
        end

        for i = 1:M
            select_val = unifrnd(0,1);
            pos = find(select_area<select_val, 1, 'last' );
            select_popV(i,:) = parent_popV(pos,:);
        end
    
        %create the new population by the selected chromosomes
        crossover_fraction = 0.8;
        crossover_count = crossover_fraction*(N-elite_count);
        mutation_count = N - elite_count - crossover_count;
    
        %a)elite children: already added
        
        %b)children created by mutation
        for i = 1:mutation_count
           children_popV(elite_count+i,:) = mutate(select_popV,curr_gen);
        end
    
        %c)children created by crossover
        for i = 1:crossover_count/2
           [children_popV(elite_count+mutation_count+2*i-1,:), children_popV(elite_count+mutation_count+2*i,:)] = crossover(select_popV); 
        end
    
        %calculate the fitness function values for the new points
        for i = 1:N
            [xV,tV,uV] = simclosedloop(children_popV(i,1),children_popV(i,2),children_popV(i,3),children_popV(i,4),children_popV(i,5),children_popV(i,6),tf);

            x1V = xV(:,1);
            x3V = xV(:,3);
            u1V = uV(:,1);
            u2V = uV(:,2);

            fit_val(i) = fitness_function(x1V,x3V,u1V,u2V,tV);
        end
        
        [elit_child, elit_child_indx] = min(fit_val);
        
        if abs(elit_child-elit_par)<1e-7
            curr_stall_gen = curr_stall_gen + 1;
        else
            curr_stall_gen = 0;
        end
    end
    
    %returned values
    chromosome = children_popV(elit_child_indx,:);
    min_fit = fit_val(elit_child_indx);
    
    [xV,tV,uV] = simclosedloop(children_popV(elit_child_indx,1),children_popV(elit_child_indx,2),children_popV(elit_child_indx,3),children_popV(elit_child_indx,4),children_popV(elit_child_indx,5),children_popV(elit_child_indx,6),tf);
    
    fprintf('\nTotal number of created generations = %d\n',curr_gen);
    fprintf('Number of stall generations = %d\n',curr_stall_gen);
    fprintf('\nMinimum value of fitness function = %f\n',min_fit);
    fprintf('\nOptimal values for the parameters of the PID controller:\n');
    fprintf('Kp1 = %f\n',chromosome(1));
    fprintf('Kp2 = %f\n',chromosome(2));
    fprintf('Ki1 = %f\n',chromosome(3));
    fprintf('Ki2 = %f\n',chromosome(4));
    fprintf('Kd1 = %f\n',chromosome(5));
    fprintf('Kd2 = %f\n',chromosome(6));
    
    %plots of x,u,e to visualise the results
    figure()
    plot(tV,xV(:,1),'b');
    yd1 = 90*pi/180 + 30*pi*cos(tV)/180;
    hold on
    plot(tV,yd1,'r');
    title('Trajectory of x1(t)')
    legend('x1(t)','yd1(t)');
    xlabel('time(s)');
    ylabel('y1(t)');
    
    figure()
    plot(tV,xV(:,3),'b');
    yd2 = 90*pi/180 + 30*pi*sin(tV)/180;
    hold on
    plot(tV,yd2,'r');
    title('Trajectory of x3(t)')
    legend('x3(t)','yd2(t)');
    xlabel('time(s)');
    ylabel('y2(t)');
    
    e1V = xV(:,1) - pi/2 - pi/6*cos(tV);
    e2V = xV(:,3) - pi/2 - pi/6*sin(tV);
    
    figure()
    plot(tV,e1V,'b');
    title('Output error e1(t)')
    xlabel('time(s)');
    ylabel('e1(t)');
    hold on
    ax = axis;
    plot([ax(1) ax(2)],(pi/180)*[1 1],'ro-');
    ylim([-0.03 0.03]);
    
    figure()
    plot(tV,e2V,'b');
    title('Output error e2(t)')
    xlabel('time(s)');
    ylabel('e2(t)');
    hold on
    ax = axis;
    plot([ax(1) ax(2)],-(pi/180)*[1 1],'ro-');
    ylim([-0.03 0.03]);
    
    figure()
    plot(tV,uV(:,1),'b');
    title('Input u1(t)')
    xlabel('time(s)');
    ylabel('u1(t)');
    hold on
    ax = axis;
    plot([ax(1) ax(2)],18*[1 1],'ro-');
    ylim([-5 20]);
    
    figure()
    plot(tV,uV(:,2),'b');
    title('Input u2(t)')
    xlabel('time(s)');
    ylabel('u2(t)');
    hold on
    ax = axis;
    plot([ax(1) ax(2)],18*[1 1],'ro-');
    ylim([-5 20]);
    
end

function f = fitness_function(x1V,x3V,u1V,u2V,t)
%returns the value of the fitness function for given x,u

    %values of J-indexes 
    J = zeros(1,6);
    
    %maximum acceptable values for each J
    Jmax = [pi/180.0 pi/180.0 pi/180.0 18 160];
    
    %matrix containing constants necessary for the fitness function
    aV = [1 1 1 200 2000];

    %calculate e
    N = length(x1V);
    
    e1V = x1V - pi/2 - pi/6*cos(t);
    e2V = x3V - pi/2 - pi/6*sin(t);
    
    %calculate J1
    if abs(e1V(N)) >= abs(e2V(N))
        J(1) = abs(e1V(N));
    else
        J(1) = abs(e2V(N));
    end
        
    %calculate J2
    max_e1 = abs(max(e1V(1:N-1)) - e1V(N));
    max_e2 = abs(e2V(N) - min(e2V(1:N-1)));
    
    if max_e1 >= max_e2
        J(2) = max_e1;
    else
        J(2) = max_e2;
    end
    
    %calculate J3
    max_abs_e1 = max(abs(e1V(10:N)));
    max_abs_e2 = max(abs(e2V(10:N)));
    
    if max_abs_e1 >= max_abs_e2
        J(3) = max_abs_e1;
    else
        J(3) = max_abs_e2;
    end
    
    %calculate J4
    max_abs_u1 = max(abs(u1V));
    max_abs_u2 = max(abs(u2V));
    
    if max_abs_u1 >= max_abs_u2
        J(4) = max_abs_u1;
    else
        J(4) = max_abs_u2;
    end
    
    %calculate J5
    J(5) = -1;
    
    for i = 2:N
        diff1 = abs((u1V(i)-u1V(i-1))/(t(i)-t(i-1)));
        diff2 = abs((u2V(i)-u2V(i-1))/(t(i)-t(i-1)));
        
        if diff1 > J(5)
            J(5) = diff1;
        elseif diff2 > J(5)
            J(5) = diff2;
        end
    end
    
    %calculate J6
    max_diff_u1 = max(u1V) - min(u1V);
    max_diff_u2 = max(u2V) - min(u2V);
    
    if max_diff_u1 >= max_diff_u2
        J(6) = max_diff_u1;
    else
        J(6) = max_diff_u2;
    end
    
    checkV = check_restrictions(J(1),J(2),J(3),J(4),J(5));
    
    %if all restrictions are satisfied
    if sum(checkV) == 0
        f = (50*J(1))/(6*Jmax(1)) + (50*J(2))/(6*Jmax(2)) + (50*J(3))/(6*Jmax(3)) + (50*J(4))/(6*Jmax(4)) + (50*J(5))/(6*Jmax(5)) + (50*J(6))/(6*(J(6)+100));    
    %if not all restrictions are satisfied
    else
        indx = find(checkV==1);
        
        f = 50;
        
        for i = 1:length(indx)
            f = f + 10*J(indx(i))/(J(indx(i))+aV(indx(i)));
        end
    end
        
end

function checkV = check_restrictions(J1,J2,J3,J4,J5)
%checks if given J satisfy the restrictions

    checkV = zeros(1,5);
    
    if J1 > pi/180.0
        checkV(1) = 1;
    end
    
    if J2 > pi/180.0
        checkV(2) = 1;
    end
    
    if J3 > pi/180.0
        checkV(3) = 1;
    end
    
    if J4 > 18
        checkV(4) = 1;
    end
    
    if J5 > 160
        checkV(5) = 1;
    end

end

function mutantV = mutate(select_popV,curr_gen)
%mutates a random chromosome from the selected population

    %M = number of chromosomes in the selected population
    %N = number of variables n each chromosome
    [M,N] = size(select_popV);
    
    mutantV = zeros(1,N);

    %mutation probability
    mut_prob = 0.01;
    
    %mutation parameters
    s = unifrnd(0,1);
    max_gens = 100*N;
    b = 6; 
    LB = 0;
    UB = [100 100 10 10 100 100];
       
    %select a random chromosome from the selected population
    rnd_indx = unidrnd(M);
    
    %make a random mutation
    for i = 1:N
        mut_rnd = unidrnd(1000);
        if mut_rnd <= 1000*mut_prob
            LBM = max([LB select_popV(rnd_indx,i)-(UB(i)-LB)*(1-s^((1-curr_gen/max_gens)^b))]);
            UBM = min([UB(i) select_popV(rnd_indx,i)+(UB(i)-LB)*(1-s^((1-curr_gen/max_gens)^b))]);
            mutantV(i) = unifrnd(LBM,UBM);
        else
            mutantV(i) = select_popV(rnd_indx,i);
        end
    end
end

function [child1,child2] = crossover(select_popV)
%creates an offspring chromosome using crossover on two parents

    %M = number of chromosomes in the selected population
    %N = number of variables n each chromosome
    [M,N] = size(select_popV);

    %select random parents (make sure they are different)
    indx1 = unidrnd(M);
    
    indx2 = indx1;
    while indx2 == indx1
        indx2 = unidrnd(M);
    end
    
    parent1 = select_popV(indx1,:);
    parent2 = select_popV(indx2,:);
    
    %choose the random cutting point
    cut_point = unidrnd(N-1);
    
    %create the two offsprings
    child1 = [parent1(1:cut_point) parent2(cut_point+1:end)];
    child2 = [parent2(1:cut_point) parent1(cut_point+1:end)];

end