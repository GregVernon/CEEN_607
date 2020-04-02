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
                        J = cell(length(P),1);
                        S = cell(length(P),1);
                        G = cell(length(P),1);
                        B = cell(length(P),1);
                        tangent = cell(length(P),1);
                        normal = cell(length(P),1);
                    end
                    if qID == 0
                        P{qp} = feElement.map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).BasisFunction{qp},[obj.Reference.Nodes.Coordinates]');
                        J{qp} = feElement.jacobian_map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).GradientBasisFunction{qp} ,[obj.Reference.Nodes.Coordinates]');
                        S{qp} = feElement.compute_integral_scaling(J{qp}, qID);
                        G{qp} = obj.Parametric.Quadrature(qID+1).GradientBasisFunction{qp} * inv(J{qp});
                        B{qp} = feElement.assemble_strain_displacement(G{qp});
                        tangent = [];
                        normal = [];
                    elseif qID == 1 || qID == 2
                        P{qp} = feElement.map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).BasisFunction{qp},[obj.Reference.Nodes.Coordinates]');
                        J{qp} = feElement.jacobian_map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).GradientBasisFunction{qp} ,[obj.Reference.Nodes.Coordinates]');
                        S{qp} = feElement.compute_integral_scaling(J{qp}, qID);
                        G{qp} = obj.Parametric.Quadrature(qID+1).GradientBasisFunction{qp} * inv(J{qp});
                        B = [];
                        tangent{qp} = feElement.extract_tangent_from_jacobian(J{qp}, qID);
                        if     qID == 1
                            normal{qp} = feElement.compute_normal([0; 0; 1], [tangent{qp}; 0]);
                        elseif qID == 2
                            normal{qp} = feElement.compute_normal([tangent{qp}; 0],[0; 0; 1]);
                        end
                        normal{qp} = normal{qp}(1:2);
                    elseif qID == 3 || qID == 4
                        P{qp} = feElement.map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).BasisFunction{qp},[obj.Reference.Nodes.Coordinates]');
                        J{qp} = feElement.jacobian_map_parametric_to_reference(obj.Parametric.Quadrature(qID+1).GradientBasisFunction{qp} ,[obj.Reference.Nodes.Coordinates]');
                        S{qp} = feElement.compute_integral_scaling(J{qp}, qID);
                        G{qp} = obj.Parametric.Quadrature(qID+1).GradientBasisFunction{qp} * inv(J{qp});
                        B = [];
                        tangent{qp} = feElement.extract_tangent_from_jacobian(J{qp}, qID);
                        if     qID == 3
                            normal{qp} = feElement.compute_normal([tangent{qp}; 0], [0; 0; 1]);
                        elseif qID == 4
                            normal{qp} = feElement.compute_normal([0; 0; 1], [tangent{qp}; 0]);
                        end
                        normal{qp} = normal{qp}(1:2);
                    end
                end
                
                obj.Reference.Quadrature(qID+1).Coordinates = P;
                obj.Reference.Quadrature(qID+1).JacobianMatrix = J;
                obj.Reference.Quadrature(qID+1).IntegralScaling = S;
                obj.Reference.Quadrature(qID+1).GradientBasisFunction = G;
                obj.Reference.Quadrature(qID+1).StrainDisplacementMatrix = B;
                obj.Reference.Quadrature(qID+1).NormalVector = normal;
            end
        end
    end
      
    methods
        function [K, KN] = compute_local_stiffness_matrix(obj)
            num_nodes = length(obj.NodeConnectivity);
            num_node_dofs = size(obj.DOFConnectivity,1);
            num_quad_points = length(obj.Parametric.Quadrature(1).Weights);
            KN = repmat({zeros(num_node_dofs,num_node_dofs)},num_nodes,num_nodes);
            D = obj.Reference.MaterialConstitutiveMatrix;
            for qp = 1:num_quad_points
                S = obj.Reference.Quadrature(1).IntegralScaling{qp};
                W = obj.Parametric.Quadrature(1).Weights{qp};
                for n1 = 1:num_nodes
                    B1 = obj.Reference.Quadrature(1).StrainDisplacementMatrix{qp}{n1};
                    for n2 = 1:num_nodes
                        B2 = obj.Reference.Quadrature(1).StrainDisplacementMatrix{qp}{n2};
                        KN{n1,n2} = KN{n1,n2} + (transpose(B1) * D * B2) * S * W;
                    end
                end
            end
            K = cell2mat(KN);
        end
        
        function [F, FN] = compute_local_internalforce_vector(obj, d_coeff)
            num_nodes = length(obj.NodeConnectivity);
            num_node_dofs = size(obj.DOFConnectivity,1);
            num_quad_points = length(obj.Parametric.Quadrature(1).Weights);
            FN = repmat({zeros(num_node_dofs,1)},num_nodes,1);
            virtual_strain = obj.compute_virtual_strain(d_coeff);
            virtual_stress = obj.compute_virtual_stress(virtual_strain);
            for qp = 1:num_quad_points
                S = obj.Reference.Quadrature(1).IntegralScaling{qp};
                W = obj.Parametric.Quadrature(1).Weights{qp};
                for n1 = 1:num_nodes
                    B = obj.Reference.Quadrature(1).StrainDisplacementMatrix{qp}{n1};
                    FN{n1} = FN{n1} + (B * virtual_stress{qp}) * S * W;
                end
            end
            F = cell2mat(FN);
        end
        
        function virtual_strain = compute_virtual_strain(obj,d_coeff)
            num_nodes = length(obj.NodeConnectivity);
            num_quad_points = length(obj.Parametric.Quadrature(1).Weights);
            virtual_strain = cell(num_quad_points,1);
            for qp = 1:num_quad_points
                for n1 = 1:num_nodes
                    B = obj.Reference.Quadrature(1).StrainDisplacementMatrix{qp}{n1};
                    if n1 == 1
                        virtual_strain{qp} = zeros(size(B,1),1);
                    end
                    virtual_strain{qp} = virtual_strain{qp} + (B * d_coeff{n1});
                end
            end
        end
        
        function virtual_stress = compute_virtual_stress(obj,virtual_strain)
            num_quad_points = length(obj.Parametric.Quadrature(1).Weights);
            virtual_stress = cell(num_quad_points,1);
            D = obj.Reference.MaterialConstitutiveMatrix;
            for qp = 1:num_quad_points
                virtual_stress{qp} = D * virtual_strain{qp};
                virtual_stress{qp} = virtual_stress{qp}(1:2);
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
        
        function refJacobian = jacobian_map_parametric_to_reference(parGradBasisFun, refNodes)
            numBasis = size(parGradBasisFun,1);
            refJacobian = zeros(size(parGradBasisFun,2), size(parGradBasisFun,2));
            for ii = 1:size(parGradBasisFun,2)
                for jj = 1:size(parGradBasisFun,2)
                    for n = 1:numBasis
                        refJacobian(ii,jj) = refJacobian(ii,jj) + parGradBasisFun(n,jj) * refNodes(n,ii);
                    end
                end
            end
        end
        
        function intScaleFactor = compute_integral_scaling(Jacobian, qID)
            if qID == 0
                intScaleFactor = det(Jacobian);
            elseif qID == 1 || qID == 2
                intScaleFactor = norm(Jacobian * [0; 1]);
            elseif qID == 3 || qID == 4
                intScaleFactor = norm(Jacobian * [1; 0]);
            end
        end
        
        function strainDispMatrix = assemble_strain_displacement(GradBasisFun)
            num_basis_fun = size(GradBasisFun,1);
            strainDispMatrix = cell(num_basis_fun,1);
            for n = 1:num_basis_fun
                strainDispMatrix{n} = [GradBasisFun(n,1), 0;
                                       0                , GradBasisFun(n,2);
                                       GradBasisFun(n,2), GradBasisFun(n,1)];
            end
        end
        
        function normal = compute_normal(V1, V2)
            normal = cross(V1, V2) / norm(cross(V1, V2));
        end
        
        function tangent = extract_tangent_from_jacobian(Jacobian, qID)
            if     qID == 1 || qID == 2
                tangent = Jacobian * [0; 1];
            elseif qID == 3 || qID == 4
                tangent = Jacobian * [1; 0];
            end
        end
        
        function D = compute_material_consitutive_matrix(E,v)
            D = (E / ((1+v)*(1-v))) * ...
                [1-v , v   , 0;
                 v   , 1-v , 0;
                 0   , 0   , (1-2*v)/2];
        end
    end
end

