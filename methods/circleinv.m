function result_full = circleinv(A,nodes,times)

result = zeros(size(A));
dim = nodes*ones(times,1);
A = mat2cell(A,dim,dim);
result = mat2cell(result,dim,dim);
    for i = 1:times
        row = times-i + 1;
        
        %compute the diagonal
        diag_result  = ones(nodes,nodes)./A{row,row};
        diag_result(isinf(diag_result)|isnan(diag_result)) = 0;
        result{row,row} = diag_result;
        
        for j = row+1:times
            temp = zeros(nodes,nodes);
            for k = row+1:j
               %[row,k,row,k,k,j]
                temp = temp - A{row,k}.*result{k,j};
            end
            
            temp= temp./A{row,row};
            temp(isinf(temp)|isnan(temp)) = 0;
            result{row,k}  = temp;
        end


    end
% A_density = cell2mat(A);
% A_density = A_density./A_density;
% A_density(isinf(A_density)|isnan(A_density)) = 0;

try
    result_full = cell2mat(result);
catch
    result_full = cell2sym(result);
end
% result = result.*A_density;
result_full = result_full;
    
end
