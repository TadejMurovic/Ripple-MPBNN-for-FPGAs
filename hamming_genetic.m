function ordered = hamming_genetic(set, iter_nb, pop_nb)
    
    % Parameters
    n0       = size(set,1);
    a        = meshgrid(1:n0);
    dist     = reshape((sum(xor(set(a,:),set(a',:)),2)),n0,n0);
    [N,dims] = size(set);
    [nr,nc]  = size(dist);
    nf = N;
    
    
    pop = zeros(pop_nb,nf);
    pop(1,:) = (1:nf);
    for k = 2:pop_nb
        pop(k,:) = randperm(nf);
    end
    
           
    global_min   = Inf;
    total_dist   = zeros(1,pop_nb);
    dist_history = zeros(1,iter_nb+1);
    tmp_pop = zeros(4,nf);
    new_pop = zeros(pop_nb,nf); 
    
    figure(1000);
    for iter = 1:iter_nb
        % Evaluate Each Population Member (Calculate Total Distance)
        for p = 1:pop_nb
            d = dist(pop(p,nf),pop(p,1)); % Closed Path
            for k = 2:nf
                d = d + dist(pop(p,k-1),pop(p,k));
            end
            total_dist(p) = d;
            if (p == 1 && iter == 1)
               dist_history(1) = d; 
            end
        end
        
        % Find the Best Route in the Population
        [min_dist,index] = min(total_dist);
        if iter == 1
            best_first = sum(dist(1,:));
        end
        dist_history(iter+1) = min_dist;
        if min_dist < global_min 
            global_min = min_dist;
            optRoute = pop(index,:);
            subplot(3,1,1);
            image(set'*100); 
            title(['Input Configuration']);
            subplot(3,1,2);
            image(set(optRoute,:)'*100); 
            title(['Best Configuration']);
            subplot(3,1,3);
            plot(1:iter+1,dist_history(1:iter+1),'r','linewidth',2);
            axis([1 iter_nb -inf best_first]);
            Improvement = round(100*(1-dist_history(iter+1)/dist_history(1)));
            title(['Best Distance: ',num2str(min_dist),', Improvement: ',num2str(Improvement),'%']);
            drawnow;
        end
                       
        % Genetic Algorithm Operators
        random_order = randperm(pop_nb);
        for p = 4:4:pop_nb
            rtes = pop(random_order(p-3:p),:);
            dists = total_dist(random_order(p-3:p));
            [ignore,idx] = min(dists); %#ok
            bestOf4Route = rtes(idx,:);
            routeInsertionPoints = sort(ceil(nf*rand(1,2)));
            I = routeInsertionPoints(1);
            J = routeInsertionPoints(2);
            for k = 1:4 % Mutate the Best to get Three New Routes
                tmp_pop(k,:) = bestOf4Route;
                switch k
                    case 2 % Flip
                        tmp_pop(k,I:J) = tmp_pop(k,J:-1:I);
                    case 3 % Swap
                        tmp_pop(k,[I J]) = tmp_pop(k,[J I]);
                    case 4 % Slide
                        tmp_pop(k,I:J) = tmp_pop(k,[I+1:J I]);
                    otherwise % Do Nothing
                end
            end
            new_pop(p-3:p,:) = tmp_pop;
        end
        pop = new_pop;
               
    end
    
    subplot(3,1,1);
    image(set'*100); 
    title(['Input Configuration']);
    subplot(3,1,2);
    image(set(optRoute,:)'*100); 
    title(['Best Configuration @ Iteration:',num2str(iter),'']);
    subplot(3,1,3);
    plot(1:iter+1,dist_history(1:iter+1),'r','linewidth',2);
    axis([1 iter_nb -inf best_first]);
    Improvement = round(100*(1-dist_history(iter+1)/dist_history(1)));
    title(['Best Distance: ',num2str(min_dist),', Improvement: ',num2str(Improvement),'%']);
    drawnow;
    
    ordered = optRoute;
    
end
