function A = dd(A,nodes,times)

dim = nodes*ones(times,1);
A = mat2cell(A,dim,dim);

for i = 1:times
    for j =1:times
        
        A{i,j} = diag(diag(A{i,j}));
        
    end
    
end

try
A = cell2mat(A);
catch
    A = cell2sym(A);
end


end
