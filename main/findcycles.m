function c = findcycles(g)
%FINDCYCLES Finds all simple cycles in graph G
%   C = FINDCYCLES(G) returns a cell array containing all simple cycles
%   found in graph G. A simple cycle is one where no node is repeated. For
%   example [1 4 7 1] is a simple cycle. However, [1 4 7 4 2 1] is a cycle,
%   but not a simple cycle.
%   C is a cell array, where each element is a vector of node indices
%   specifying a cycle in G. Each path starts with its lowest-index node.
%
%   Note I just quickly coded this up to answer a question, no guarantees.

c = {};

for ii=1:numnodes(g)
    % Find all cycles starting with node ii, which only contain nodes
    % with indices >= ii.
    c = findcycleRecursive(ii, g, c, ii);
end

end

function c = findcycleRecursive(currentPath, g, c, minNodeID)

if isa(g, 'graph')
    successorLast = reshape(neighbors(g, currentPath(end)), 1, []);
else
    successorLast = reshape(successors(g, currentPath(end)), 1, []);
end

for s = successorLast
    if s < minNodeID
        % This has already been found with starting point < minNodeID
        % (avoid returning [2 4 1], [1 2 4], [4 1 2] -> always start with
        % smallest number)
        continue;
    end
    
    if s == currentPath(1)
        % Found a cycle
        c{end+1} = [currentPath s]; %#ok<AGROW>
    elseif any(currentPath == s)
        % Contains a subcycle, this would result in a cycle that is not
        % simple at best.
    else
        % Add node s and continue growing the path.
        c = findcycleRecursive([currentPath s], g, c, minNodeID);
    end
end

end