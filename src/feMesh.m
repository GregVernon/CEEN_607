classdef feMesh
    %feMesh Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NodeConnectivity
        DOFConnectivity
        ElementSets feElementSet
        SurfaceSets feSurfaceSet
        NodeSets    feNodeSet
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
                K(local_2_global_dof,local_2_global_dof) = K(local_2_global_dof,local_2_global_dof) + k;
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
                Fext(local_2_global_dof) = Fext(local_2_global_dof) + f;
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
                Fint(local_2_global_dof) = Fint(local_2_global_dof) + f;
            end
        end
        
        function obj = assignBoundaryConditions(obj, BC)
            num_boundary_condtions = length(BC);
            for bndy_id = 1:num_boundary_condtions
                if strcmpi(BC(bndy_id).Type, "dirichlet")
                    mesh_node_set_id = strcmpi(BC(bndy_id).NodeSet.Name, {obj.NodeSets(:).Name});
                    NodeSet = obj.NodeSets(mesh_node_set_id);
                    num_nodes_in_NodeSet = length(NodeSet.GlobalNodeID);
                    for n = 1:num_nodes_in_NodeSet
                        parent_elems = NodeSet.ElementID{n};
                        for e = 1:length(parent_elems)
                            global_elem_id = parent_elems(e);
                            loc_node_id = find(obj.Elements(global_elem_id).NodeConnectivity == NodeSet.GlobalNodeID(n));
                            num_local_dof = size(obj.Elements(global_elem_id).DOFConnectivity,1);
                            loc_dof_id = num_local_dof * (loc_node_id - 1) + BC(bndy_id).DOF;
                            obj.Elements(global_elem_id).Reference.DirichletConditions{loc_node_id}(BC(bndy_id).DOF) = BC(bndy_id).Value;
                        end
                    end
                end
            end
        end
        
        function obj = assignLoadConditions(obj, LC)
            num_load_conditions = length(LC);
            for load_id = 1:num_load_conditions
                if     strcmpi(LC(load_id).Type, "body")
                    mesh_elem_set_id = strcmpi(LC(load_id).ElementSet.Name, {obj.ElementSets(:).Name});
                    ElementSet = obj.SurfaceSets(mesh_elem_set_id);
                    ElemID = ElementSet.ElementID;
                    M = LC(load_id).Magnitude;
                    D = LC(load_id).Direction;
                    BF = M * (D ./ norm(D));
                    for loc_elem_id = 1:length(ElemID)
                        global_elem_id = ElemID(loc_elem_id);
                        obj.Elements(global_elem_id).Reference.BodyForce = BF;
                    end
                elseif strcmpi(LC(load_id).Type, "pressure")
                    mesh_surf_set_id = strcmpi(LC(load_id).SurfaceSet.Name, {obj.SurfaceSets(:).Name});
                    SurfaceSet = obj.SurfaceSets(mesh_surf_set_id);
                    ElemID = SurfaceSet.ElementID;
                    M = LC(load_id).Magnitude;
                    for loc_elem_id = 1:length(ElemID)
                        % Get global element id
                        global_elem_id = ElemID(loc_elem_id);
                        % Get local element side id
                        loc_side_id = SurfaceSet.LocalSideID(loc_elem_id);
                        % Assign pressure to local element side id
                        obj.Elements(global_elem_id).Reference.SurfacePressure{loc_side_id} = M;
                    end
                elseif strcmpi(LC(load_id).Type, "force")
                    if     isempty(LC(load_id).SurfaceSetName) == false
                        % Node Force defined on all nodes in a surface set
                        mesh_surf_set_id = strcmpi(LC(load_id).SurfaceSet.Name, {obj.SurfaceSets(:).Name});
                        SurfaceSet = obj.SurfaceSets(mesh_surf_set_id);
                        SurfaceSet = LC(load_id).SurfaceSet;
                        ElemID = SurfaceSet.ElementID;
                        M = LC(load_id).Magnitude;
                        D = LC(load_id).Direction;
                        BF = M * (D ./ norm(D));
                        for loc_elem_id = 1:length(ElemID)
                            global_elem_id = ElemID(loc_elem_id);
                            loc_node_id = NodeSet.LocalNodeID(loc_elem_id);
                            obj.Elements(global_elem_id).Reference.NodeForce{loc_node_id} = BF;
                        end
                    elseif isempty(LC(load_id).NodeSetName) == false
                        % Node Force defined on all nodes in a nosde set
                        mesh_node_set_id = strcmpi(LC(load_id).NodeSet.Name, {obj.NodeSets(:).Name});
                        NodeSet = obj.SurfaceSets(mesh_node_set_id);
                        ElemID = NodeSet.ElementID;
                        M = LC(load_id).Magnitude;
                        D = LC(load_id).Direction;
                        BF = M * (D ./ norm(D));
                        for loc_elem_id = 1:length(ElemID)
                            global_elem_id = ElemID(loc_elem_id);
                            loc_node_id = NodeSet.LocalNodeID(loc_elem_id);
                            obj.Elements(global_elem_id).Reference.NodeForce{loc_node_id} = BF;
                        end
                    end
                end
            end
        end
    end
    
    methods
        function plotMesh(obj)
            hold on;
            num_elems = length(obj.Elements);
            for e = 1:num_elems
                xy = [obj.Elements(e).Reference.Nodes.Coordinates];
                x = xy(1,[1 2 4 3]);
                y = xy(2,[1 2 4 3]);
                patch('XData',x,'YData',y,'FaceColor','None','EdgeColor','k','LineWidth',1,'Marker','.','MarkerSize',20)
            end
            
            for e = 1:num_elems
                xy = [obj.Elements(e).Reference.Nodes.Coordinates];
                num_elem_nodes = length(x);
                for n = 1:num_elem_nodes
                    global_node_id = obj.Elements(e).NodeConnectivity(n);
                    text(xy(1,n),xy(2,n),num2str(global_node_id),'FontSize',24)
                end
            end
        end
    end
end

