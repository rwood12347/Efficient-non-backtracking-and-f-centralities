function d = Centrality_Measure(As,M, measure, alpha)

% Centrality_Measure gives the centrality of
% each node in a graph with the adjacency matrix A. 

%DONE: check whether or not we're inputting a single matrix or a list of
%matrices

%TODO: Implement NBT.


no_t = size(As,2);
nodes = size(As{end},2);
I = eye(size(M));


% Options for the measure:
%
% 'deg'         : Degree centrality
% 'Katz_half'   : Katz parameter used to study similarity in texts
% 'Katz_Google' : Katz parameter used in Google's page-rank 
% 'Katz_deg'    : alpha = 1/(||A||_inf+1)
% 'Katz_min'    : alpha = (1-exp(-lamda))/lamda
% 'exp'         : Total subgrap communicability
% 'Sub_cent'    : Subgraph centrality
% 'between'     : Betweennes centrality
% 'NBTW'        : Nonbacktracking walk



switch measure
    
    case 'deg'
        %Degree centrality for a single time frame.
        d = zeros(nodes,1);
            for i = 1:size(As,2)
                d = d + As{i}*ones(nodes,1);
            end
            

		 case 'Katz_deg'
        alpha =  1/max(cellfun(@(x)norm(x,'inf')+1,As));
            d = ones(nodes,1);
            for i = 1:size(As,2)
                d = (eye(nodes) - alpha*As{no_t - i + 1}) \ d;
            end
        
        
    case 'Katz_min'
        lambda = cellfun(@(x) eigs(x,1), As);
        lambda = max(abs(lambda));
        alpha = (1-exp(-lambda))/lambda;
        if isnan(alpha)
            alpha = 1;
        end
            d = ones(nodes,1);
            for i = 1:size(As,2)
                d = (eye(nodes) - alpha*As{no_t - i + 1}) \ d;
            end
            %sometimes matlab gets confused and adds minuses
            d = abs(d);
    case 'Katz'
            lambda = cellfun(@(x) eigs(x,1), As);
            lambda = max(abs(lambda));
            alpha=alpha/lambda;
            d = ones(nodes,1);
            for i = 1:size(As,2)
                d = (eye(nodes) - alpha*As{no_t - i + 1}) \ d;
            end
    case 'KatzWindow'
            lambda = max(abs(eigs(As{end},1)));
            alpha = alpha/lambda;
            d = ones(nodes,1);
            d = (eye(nodes) - alpha*As{end}) \ d;
    case 'NBTWWindow'
            lambda = max(abs(eigs(As{end},1)));
            alpha = alpha/lambda;
            d  = NewnodeNBTW(As{end}, alpha, nodes, 1);
            d = d*ones(nodes,1);
            
    case 'exp'
%         [nodes, no_t]
%         [T,B] = balance(full(M*alpha));
%         test = T*expm(B)/T;
%         test = test*ones(size(test,2),1);
        d = expm(full(M*(alpha)))*ones(size(M,1),1);

        
        d =  d(1:nodes);
        d = d/norm(d,2);
                        d = d.*(d >eps);
    case 'Sub_cent_exp'
%         tic

        d = expm(full(M*alpha))*repmat(eye(nodes), no_t, 1);
            
        d = diag(d);
                d = d/norm(d,2);
                d = d.*(d >eps);
%        oldt = toc 
%        tic
%             test = expm(full(M)*alpha);
%             test = test(1:nodes, :);
%             cent = zeros(nodes,1);
%             for i = 1:no_t
%                cent = cent + diag(test(:, (i-1)*nodes+1:i*nodes)); 
%             end
%             test = cent;
%             newt = toc
%             sum(d - test)

    case 'Sub_cent_Katz'
            lambda = cellfun(@(x) eigs(x,1), As);
            lambda = max(abs(lambda));
            alpha=alpha/lambda;
            d = eye(nodes);
            for i = 1:size(As,2)
                d = (eye(nodes) - alpha*As{no_t - i + 1})  \ d;
            end
             d = diag(d);
    case 'Sub_cent_NBTW'
        lambda=1/abs(eigs(full(M),1));
        if isinf(lambda)
            lambda = 1;
        end
         alpha=alpha*lambda;
         d = NewnodeNBTW(full(M), alpha, nodes, no_t);
%          d = d*ones(nodes*no_t,nodes);
         d = dd(d,nodes,no_t);
         d  = d(1:nodes,:)*ones(nodes*no_t,1);
    case 'NBTW'

          lambda=1/abs(eigs(full(M),1));
          if isinf(lambda)
            lambda = 1;
        end
         alpha=alpha*lambda;
         
        d = NewnodeNBTW(full(M), alpha, nodes, no_t);
        d = d(1:nodes,:)*ones(nodes*no_t,1);
       

	
        
         
        
end




%         
%         


%         
% 	case 'Sub_cent_NBTW_exp'
%         if ~isCellArray
%             d = diag(NBTW_exp(As,alpha));
%         else
%             error('Sub_cent_NBTW_exp not implemented yet')
%         end
%   	case 'NBTW_exp'
%         if ~isCellArray
%             d = NBTW_exp(As,alpha)*I;
%         else
%             error('NBTW_exp not implemented yet')
%         end
%         
%             
% 	case 'evec'
%         if ~isCellArray
%             lambda = max(eig(As));
%             d=abs(null(As-lambda*I_N));
%         else
%             error('evec not implemented yet')
%         end
% 		
%  	case 'NBTWevec'
%         if ~isCellArray
%             [~, D]=NBTW(As,42);
%             lambda=min(polyeig(I_N,-As,D-I_N));
%             d=abs(null(I_N-As*lambda+(D-I_N)*lambda^2));
%         else
%             error('NBTWevec not implemented yet')
%         end
% 	case 'between'
%         if ~isCellArray
%             d = betweenness(As);
%         else
%             error('between not implemented yet')
%         end
% end
% 
% 


