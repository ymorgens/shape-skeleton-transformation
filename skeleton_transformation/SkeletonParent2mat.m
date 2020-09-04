function mat = SkeletonParent2mat(skeleton)

nindex = length(skeleton);
mat = zeros(1,nindex);
for i = 1:nindex
mat(i) = skeleton(i).parent;
    

end