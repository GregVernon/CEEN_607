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
        function K = assemble_global_stiffness_matrix(obj)
            num_elems = length(obj.Elements);
            if isempty(obj.DOFConnectivity)
                dof_connect = obj.gather_dof_connectivity;
                num_dof = max(dof_connect(:));
            else
                num_dof = max(obj.DOFConnectivity(:));
            end
            
            K = zeros(num_dof, num_dof);
            for e = 1:num_elems
                k = obj.Elements(e).Reference.StiffnessMatrix;
                local_2_global_dof = obj.Elements(e).DOFConnectivity(:);
                K(local_2_global_dof,local_2_global_dof) = k;
            end
        end
        
        function Fext = assemble_global_externalforce_vector(obj)
            num_elems = length(obj.Elements);
            if isempty(obj.DOFConnectivity)
                dof_connect = obj.gather_dof_connectivity;
                num_dof = max(dof_connect(:));
            else
                num_dof = max(obj.DOFConnectivity(:));
            end
            
            Fext = zeros(num_dof,1);
            for e = 1:num_elems
                f = obj.Elements(e).Reference.ExternalForceVector;
                local_2_global_dof = obj.Elements(e).DOFConnectivity(:);
                Fext(local_2_global_dof) = f;
            end
        end
        
        function Fint = assemble_global_internalforce_vector(obj)
            num_elems = length(obj.Elements);
            if isempty(obj.DOFConnectivity)
                dof_connect = obj.gather_dof_connectivity;
                num_dof = max(dof_connect(:));
            else
                num_dof = max(obj.DOFConnectivity(:));
            end
            
            Fint = zeros(num_dof,1);
            for e = 1:num_elems
                f = obj.Elements(e).Reference.InternalForceVector;
                local_2_global_dof = obj.Elements(e).DOFConnectivity(:);
                Fint(local_2_global_dof) = f;
            end
        end
    end
end

