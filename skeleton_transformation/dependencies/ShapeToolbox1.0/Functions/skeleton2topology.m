function topology = skeleton2topology(skeleton,starting_index)
% topology = skeleton2topology(skeleton,starting_index)
%
% Convert an ordinary skeleton into an pure topological graph, indicating
% who is parent of whom. A set of pairs [[parent,child],...].
% Starts from starting_index (for recursive calls).
%
this_parent_children = [];
for i=1:length(skeleton)
    if skeleton(i).parent == starting_index
        this_parent_children = [this_parent_children; i];
    end;
end;
if isempty(this_parent_children)  % Bottom out when there are no children
    topology = []; 
else
    pairs_with_immediate_children = [];
    childrens_topologies = [];
    for i=1:length(this_parent_children)
        this_topology = skeleton2topology(skeleton,this_parent_children(i));
        childrens_topologies = [childrens_topologies; this_topology];
        pairs_with_immediate_children = [pairs_with_immediate_children; starting_index this_parent_children(i)];
    end;
    topology = [pairs_with_immediate_children; childrens_topologies];    
end;
if starting_index==1
    topology = [topology; -1 1]; % add root parent = -1 by convention
end;
