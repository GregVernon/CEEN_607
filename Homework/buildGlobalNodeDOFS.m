function GNode2GDOF = buildGlobalNodeDOFS(numGlobalNodes,nDOFS)
    GNode2GDOF = zeros(nDOFS,numGlobalNodes);    
    gdof = 0;
    for n = 1:numGlobalNodes
        for ldof = 1:nDOFS
            gdof = gdof + 1;
            GNode2GDOF(ldof,n) = gdof;
        end
    end
end