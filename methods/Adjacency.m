function A = Adjacency(C,adj,threshold)

% ADJACENCY constructs the adjacency matrix 
% related to the correlation matrix C.
%
% C         : Correlation matrix
% adj       : Option for the adjoint matrix.
% threshold : Threshold 
%             (An edge between two nodes exist if and only if
%             the cross-correlation between two nodes are
%             higher then a specied threshold value.)
 

n = length(C);
I = eye(n);

if nargin < 3, threshold = 0; end

switch adj
    
    case '1'
        A = C;
    case '2'
        A = C-I;
    case '3'
        A = C > threshold;
        A = double(A);
    case '4'
        A = (C-I)>threshold;
        A = double(A);
    case '5'
      A = C.*(C>threshold);

    case '6'
        A = (C-I).*((C-I)>threshold);

end
%We might end up with negative entries just because of case 2
A = abs(A);


end