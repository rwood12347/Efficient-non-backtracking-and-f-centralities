function result =  trans(A,nodes,times)



dim = nodes*ones(times,1);
%Inefficient, you could just split the blocks and since we know the
%structure of M, transpose A1...A2 and rebuild; at the moment O(t^2n^2)

for i = 1:times
    for j =1:times

        A((i-1)*(nodes)+1:i*nodes,(j-1)*nodes+1:j*nodes) =A((i-1)*nodes+1:i*nodes,(j-1)*nodes+1:j*nodes).';

    end
    
end
result = A;
end