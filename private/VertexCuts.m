function cuts = VertexCuts(E)
% assert: sparse, symmetric graph

bc = biconnected_components_mex(double(E));

cuts = false(1,length(E));
cuts(bc(bc > 0)) = true;

end