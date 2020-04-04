classdef feElementParametric
    %feElementConfiguration Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Descriptors
        Degree
        Dimension
        GlobalID
        Type
                
        % Children
        Nodes
        Quadrature
        
        % Physics
        MaterialConstitutiveMatrix
        
        % Linear Algebra Stuff
        StiffnessMatrix
        InternalForceVector
        ExternalForceVector
    end
    
    methods
        % Constructors
        function obj = feElementParametric()
            %feElementConfiguration Construct an instance of this class
        end
    end
    
    methods (Static)
        % Nodal Information
        function nodalCoords = computeNodeCoordinates(ElementDegree, loc_nod_id)
            if     length(ElementDegree) == 1
                x = linspace(-1,1,ElementDegree(1));
                nodalCoords = [x(loc_nod_id)];
            elseif length(ElementDegree) == 2
                nID = 0;
                x = linspace(-1,1,ElementDegree(1)+1);
                y = linspace(-1,1,ElementDegree(2)+1);
                for ny = 1:ElementDegree(2)+1
                    for nx = 1:ElementDegree(1)+1
                        nID = nID + 1;
                        if nID == loc_nod_id
                            nodalCoords = [x(nx); y(ny)];
                            return
                        end
                    end
                end
            elseif length(ElementDegree) == 3
                nID = 0;
                x = linspace(-1,1,ElementDegree(1));
                y = linspace(-1,1,ElementDegree(2));
                z = linspace(-1,1,ElementDegree(3));
                for nz = 1:ElementDegree(3)+1
                    for ny = 1:ElementDegree(2)+1
                        for nx = 1:ElementDegree(1)+1
                            nID = nID + 1;
                            if nID == loc_nod_id
                                nodalCoords = [x(nx); y(ny); z(nz)];
                                return
                            end
                        end
                    end
                end
            end
        end
    end
    
     methods (Static)
        % Basis Functions
        function N = LagrangeBasis(xi, Degree)
            nPts = Degree+1;
            nDim = length(Degree);
            if    nDim == 1
                if Degree == 1
                    N = [(1-xi)/2;
                        (1+xi)/2];
                elseif Degree == 2
                    N = [xi * (xi/2 - 1/2);
                        -(xi - 1) * (xi + 1);
                        xi * (xi/2 + 1/2)];
                end
            elseif nDim == 2
                N = zeros(prod(nPts),1);
                N1 = feElementParametric.LagrangeBasis(xi(1), Degree(1));
                N2 = feElementParametric.LagrangeBasis(xi(2), Degree(2));
                n = 0;
                for jj = 1 : nPts(2)
                    for ii = 1 : nPts(1)
                        n = n+1;
                        N(n) = N1(ii) * N2(jj);
                    end
                end
            elseif nDim == 3
                N = zeros(prod(nPts));
                N1 = feElementParametric.LagrangeBasis(xi(1), Degree(1));
                N2 = feElementParametric.LagrangeBasis(xi(2), Degree(2));
                N3 = feElementParametric.LagrangeBasis(xi(3), Degree(3));
                n = 0;
                for kk = 1 : nPts(3)+1
                    for jj = 1 : nPts(2)+1
                        for ii = 1 : nPts(1)+1
                            n = n+1;
                            N(n) = N1(ii) * N2(jj) * N3(kk);
                        end
                    end
                end
            end
        end
        
        function DN = Grad_LagrangeBasis(xi,Degree)
            nDim = length(Degree);
            if    nDim == 1
                if     Degree == 1
                    DN = [-1/2
                        +1/2];
                elseif Degree == 2
                    DN = [xi - 1/2;
                        -2 * xi;
                        xi + 1/2];
                end
            elseif nDim == 2
                if     Degree == 1
                    DN = [+xi(2)/4 - 1/4, +xi(1)/4 - 1/4;
                          -xi(2)/4 + 1/4, -xi(1)/4 - 1/4;
                          -xi(2)/4 - 1/4, -xi(1)/4 + 1/4;
                          +xi(2)/4 + 1/4, +xi(1)/4 + 1/4];
                elseif Degree == 2
                    DN = [(xi(2)*(2*xi(1) - 1)*(xi(2) - 1))/4  , (xi(1)*(2*xi(2) - 1)*(xi(1) - 1))/4;
                          -xi(1)*xi(2)*(xi(2) - 1)             , -((xi(1)^2 - 1)*(2*xi(2) - 1))/2;
                          (xi(2)*(2*xi(1) + 1)*(xi(2) - 1))/4  , (xi(1)*(2*xi(2) - 1)*(xi(1) + 1))/4;
                          -((2*xi(1) - 1)*(xi(2)^2 - 1))/2     , -xi(1)*xi(2)*(xi(1) - 1);
                          2*xi(1)*(xi(2)^2 - 1)                , 2*xi(2)*(xi(1)^2 - 1);
                          -((2*xi(1) + 1)*(xi(2)^2 - 1))/2     , -xi(1)*xi(2)*(xi(1) + 1);
                          (xi(2)*(2*xi(1) - 1)*(xi(2) + 1))/4  , (xi(1)*(2*xi(2) + 1)*(xi(1) - 1))/4;
                          -xi(1)*xi(2)*(xi(2) + 1)             , -((xi(1)^2 - 1)*(2*xi(2) + 1))/2;
                          (xi(2)*(2*xi(1) + 1)*(xi(2) + 1))/4  , (xi(1)*(2*xi(2) + 1)*(xi(1) + 1))/4];
                end
            elseif nDim == 3
                if     Degree == 1
                elseif Degree == 2
                end
            end
        end
    end
end

