function M = buildM(oldM, newA)
%Returns the new matrix M = [A1 A2; 0 A2] etc.

if nargin < 2, M = oldM, newA = oldM; end

if nargin > 1
% M = [oldM repmat(newA, size(oldM,1)/size(newA, 1),1); zeros(size(newA,1),size(oldM,2)) newA];

%test----
newsize = size(oldM) + size(newA);
M = sparse(newsize(1),newsize(2));
M(1:size(oldM,1), 1:size(oldM,2)) = oldM;
try
M(1:size(oldM,1)+ size(newA,1), size(oldM,2)+1:size(oldM,2) + size(newA,2) ) = sparse(repmat(newA, size(M,1)/size(newA, 1),1));
catch
    M
end
end

   
    


end