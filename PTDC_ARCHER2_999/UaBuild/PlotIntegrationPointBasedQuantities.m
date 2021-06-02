function [FigHandle,ColorbarHandel]=PlotIntegrationPointBasedQuantities(CtrlVar,IntTriangulation,IntCoordinates,Value,varargin)
    
    %  Produces (scalar) color plot of a quantity defined at integration points
    %
    %  To plot exx for example do:
    %  figure ; PlotIntegrationPointBasedQuantities(CtrlVar,DTintTriInside,DTint.X,exx(:));
    %
    %  if IntTriangulation and IntCoordinates are not known use:
    %    [DTxy,tri,DTint,DTintTriInside,Xint,Yint,xint,yint,Iint]=...
    %           TriangulationNodesIntegrationPoints(coordinates,connectivity,Boundary.EdgeCornerNodes,nip)
    %  and then do:
    %    figure ; PlotIntegrationPointBasedQuantities(CtrlVar,DTintTriInside,DTint.X,exx(Iint));
    %
    
    Value=Value(:);
    
    FigHandle=patch('faces',IntTriangulation,'vertices',IntCoordinates/CtrlVar.PlotXYscale,...
        'FaceVertexCData',Value,'CDataMapping','scaled','EdgeColor','none','FaceColor','interp',varargin{:}) ;
    
    ColorbarHandel=colorbar;
    axis equal tight
    
    return
    
    
end