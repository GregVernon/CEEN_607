classdef feElement
    %feElement A single finite element object
    %   Detailed explanation goes here
    
    properties
        % Topology
        NodeConnectivity
        DOFConnectivity
        
        % Descriptors
        Degree
        Dimension
        GlobalID
        Type
                
        % Configurations
        Parametric   feElementParametric
        Reference    feElementReference
        Deformed     feElementDeformed
    end
    
    methods
        % Constructors
        function obj = feElement()
            %feElement Construct an empty instance of this class
            %   Detailed explanation goes here
        end
    end
    
    methods (Static)
        % Basis Functions
        function N = LagrangeBasis(xi, Degree)
            nPts = Degree+1;
            nDim = length(Degree);
            switch nDim
                case 1
                    if Degree == 1
                        N = [(1-xi)/2;
                             (1+xi)/2];
                    elseif Degree == 2
                        N = [xi * (xi/2 - 1/2); 
                             -(xi - 1) * (xi + 1); 
                             xi * (xi/2 + 1/2)];
                    end
                case 2
                    N = zeros(prod(nPts),1);
                    N1 = feElement.LagrangeBasis(xi(1), Degree(1));
                    N2 = feElement.LagrangeBasis(xi(2), Degree(2));
                    n = 0;
                    for jj = 1 : nPts(2)
                        for ii = 1 : nPts(1)
                            n = n+1;
                            N(n) = N1(ii) * N2(jj);
                        end
                    end
                case 3
                    N = zeros(prod(nPts));
                    N1 = feElement.LagrangeBasis(xi(1), Degree(1));
                    N2 = feElement.LagrangeBasis(xi(2), Degree(2));
                    N3 = feElement.LagrangeBasis(xi(3), Degree(3));
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
    end
    
    methods
        % Change of coordinates / configuration mappings
    end
end

