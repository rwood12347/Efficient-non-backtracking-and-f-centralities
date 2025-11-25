function result = newcircle(A,B,nodes,times)

%A, B things to be circled
%n -no nodes
%t - no-times frames

%Assume A, B are both upper diagonal.

dim = nodes*ones(times,1);
result = zeros(size(A));


for i = 1:times
    
    for j = i:times
        temp_sum = zeros(nodes,nodes);
        for k = 1:j
        try    
        temp_sum = temp_sum + A((i-1)*nodes+1:i*nodes,(k-1)*nodes+1:k*nodes).*B((k-1)*nodes+1:k*nodes,(j-1)*nodes+1:j*nodes);
        catch 
            disp('wtf')
        end
        end
        try
         result((i-1)*nodes+1:i*nodes,(j-1)*nodes+1:j*nodes)= temp_sum;
        catch
            result = sym(result);
            result((i-1)*nodes+1:i*nodes,(j-1)*nodes+1:j*nodes)= temp_sum;
        end
     end
end

% try

%  result = sparse(cell2mat(result));
 
% catch
%     try
%         result = cell2sym(result);
%     catch
%     disp('(RYAN) Error in conversion to matrix');
%     end

end
