classdef feMesh
    %feMesh Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NodeConnectivity
        DOFConnectivity
        Elements    feElement
    end
    
    methods
        % Constructor
        function obj = feMesh(nElems)
            %feMesh Construct an instance of this class
            %   Detailed explanation goes here
            obj.Elements = repmat(feElement(), nElems, 1);
        end
    end
    
    methods 
        % Precomputations
        function obj = precomputeMesh(obj)
            num_elems = length(obj.Elements);
            for e = 1:num_elems
                obj.Elements(e) = obj.Elements(e).precompute();
            end
        end
    end
    
    methods
        function NodeConnect = gather_node_connectivity(obj)
            NodeConnect = [obj.Elements.NodeConnectivity];
        end
        
        function DOFConnect = gather_dof_connectivity(obj)
            DOFConnect = [obj.Elements.DOFConnectivity];
        end
    end
    
    methods
        function [K] = assemble_global_stiffness_matrix(obj)
            num_elems = length(obj.Elements);
            num_dof = max(obj.DOFConnectivity(:));
            K = zeros(num_dof, num_dof);
            for e = 1:num_elems
                k = obj.Elements(e).Reference.StiffnessMatrix;
                local_2_global_dof = obj.Elements(e).DOFConnectivity;
                K(local_2_global_dof,local_2_global_dof) = k;
            end
        end
    end
end

