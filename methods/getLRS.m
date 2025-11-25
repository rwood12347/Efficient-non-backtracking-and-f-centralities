function [L,R,B] = getLRB(A)
%-------------------------------------------------------------------------------
%Computes the source, target and weight matrix thereby allowing for
%function based NBT centralities for both static and time-evolving
%networks.

%-------------------------------------------------------------------------------
% INPUTS
%A : The weighted adjacency matrix for which we wish to return the source,
%   target and weight matrices.

%-------------------------------------------------------------------------------
% OUTPUTS
% L - Source matrix
% R - Target matrix
% S - Weight matrix
%-------------------------------------------------------------------------------


g = digraph(A);

edges = g.Edges;
no_e = size(edges,1);
nodes = g.Nodes;
no_n = size(nodes,1);
B = zeros(no_e, no_e);
L = zeros(no_e, no_n);
R = zeros(no_e, no_n);

for edge = 1:size(edges, 1)
 i = edges(edge,:).EndNodes(1);
 j = edges(edge,:).EndNodes(2);
 L(edge, i) = 1;
 R(edge, j) = 1;
 for second_edge = 1:size(edges,1)
        r = edges(second_edge,:).EndNodes(1);
        s = edges(second_edge,:).EndNodes(2);
        if and(j == r, i ~= s)
            B(edge,second_edge) = 1
        end
 end



end
end

