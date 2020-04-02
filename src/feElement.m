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
    
    methods
        % Precomputations
        function obj = precompute(obj)
            % Precompute Parametric Configuration
            num_dim = obj.Dimension;
            num_quadrature = length(obj.Parametric.Quadrature);
            for qID = 0:(num_quadrature-1)
                num_qp = obj.computeNumberOfQuadraturePoints(qID);
                [P,W] = feQuadrature.assembleQuadratureRule(num_dim, num_qp, qID);
                for qp = 1:length(P)
                    if qp == 1
                        bFun = cell(length(P),1);
                        grad_bFun = cell(length(P),1);
                    end
                    if     qID == 0
                        bFun{qp} = feElementParametric.LagrangeBasis(P{qp},obj.Degree);
                        grad_bFun{qp} = feElementParametric.Grad_LagrangeBasis(P{qp}, obj.Degree);
                    elseif qID == 1 || qID == 2
                        bFun{qp} = feElementParametric.LagrangeBasis(P{qp},obj.Degree);
                        grad_bFun{qp} = feElementParametric.Grad_LagrangeBasis(P{qp}, obj.Degree);
                    elseif qID == 3 || qID == 4
                        bFun{qp} = feElementParametric.LagrangeBasis(P{qp},obj.Degree);
                        grad_bFun{qp} = feElementParametric.Grad_LagrangeBasis(P{qp}, obj.Degree);
                    end
                end
                obj.Parametric.Quadrature(qID+1).Coordinates = P;
                obj.Parametric.Quadrature(qID+1).Weights     = W;
                obj.Parametric.Quadrature(qID+1).BasisFunction = bFun;
                obj.Parametric.Quadrature(qID+1).GradientBasisFunction = grad_bFun;
            end
            
            % Precompute Reference Configuration
            for qID = 0:(num_quadrature-1)
                num_qp = obj.computeNumberOfQuadraturePoints(qID);
                [P,W] = feQuadrature.assembleQuadratureRule(num_dim, num_qp, qID);
                for qp = 1:length(P)
                    if qp == 1
                        P = cell(length(P),1);
                    end
                    if qID == 0
                        P{qp} = feElement.map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).BasisFunction{qp},[obj.Reference.Nodes.Coordinates]');
                    elseif qID == 1 || qID == 2
                        P{qp} = feElement.map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).BasisFunction{qp},[obj.Reference.Nodes.Coordinates]');
                    elseif qID == 3 || qID == 4
                        P{qp} = feElement.map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).BasisFunction{qp},[obj.Reference.Nodes.Coordinates]');
                    end
                end
                obj.Reference.Quadrature(qID+1).Coordinates = P;
            end
        end
    end
      
    methods
        function num_qp = computeNumberOfQuadraturePoints(obj, qID)
            nDim = obj.Dimension;
            if     nDim == 1
                num_qp = ceil((obj.Degree + 1) / 2) + 1;
                
            elseif nDim == 2
                if     qID == 0  % Solid Quadrature
                    num_qp = ceil((obj.Degree + 1) ./ 2) + 1;
                elseif qID == 1 || qID == 2 % XMin and XMax Sides
                    num_qp = ceil((obj.Degree(2) + 1) / 2) + 1;
                elseif qID == 3 || qID == 4 % YMin and YMax Sides
                    num_qp = ceil((obj.Degree(1) + 1) / 2) + 1;
                end
                
            elseif nDim == 3
                if     qID == 0  % Solid Quadrature
                    num_qp = ceil((obj.Degree + 1) ./ 2) + 1;
                elseif qID == 1 || qID == 2 % XMin and XMax Sides
                    num_qp = ceil((obj.Degree([3 2]) + 1) ./ 2) + 1;
                elseif qID == 3 || qID == 4 % YMin and YMax Sides
                    num_qp = ceil((obj.Degree([1 3]) + 1) ./ 2) + 1;
                elseif qID == 3 || qID == 4 % ZMin and ZMax Sides
                    num_qp = ceil((obj.Degree([1 2]) + 1) ./ 2) + 1;
                end
            end
        end
    end
    
    methods (Static)
        function refConfig = map_parametric_to_reference(parBasisFun, refNodes)
            numBasis = size(parBasisFun,1);
            refConfig = zeros(size(refNodes,2),1);
            for n = 1:numBasis
                refConfig = refConfig + parBasisFun(n) * refNodes(n,:)';
            end
        end
    end
end

