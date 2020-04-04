classdef feSurfaceSet
    %feSurfaceSet Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        ID
        ElementID
        LocalSideID
        LocalNodeID
        GlobalNodeID
    end
    
    methods
        % Constructors
        function obj = feSurfaceSet()
            %feSurfaceSet Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = build_from_exodus(obj, EXO, setName)
            % Exodus.side_sets{:,set_id}
                % 1 - name
                % 2 - id
                % 3 - element list
                % 4 - side list
                % 5 - nodes per side list
                % 6 - node list
                % 7 - distribution factors
            set_id = find(strcmpi(setName, EXO.side_sets(1,:)));
            obj.Name = EXO.side_sets{1,set_id};
            obj.ID = EXO.side_sets{2,set_id};
            obj.ElementID = EXO.side_sets{3,set_id};
            obj.LocalSideID = feSurfaceSet.ExodusSideID_to_TensorSideID("QUAD4", EXO.side_sets{4,set_id});
            obj.GlobalNodeID = feSurfaceSet.EXO.side_sets{6,set_id};
        end
    end
    
    methods (Static)
        function TP_side_id = ExodusSideID_to_TensorSideID(elementType, EXO_sideID)
            if     strcmpi(elementType, "QUAD4")
                EXO_2_TP_side_map = [4 2 1 3];
                TP_side_id = EXO_2_TP_side_map(EXO_sideID);
            elseif strcmpi(elementType, "HEX8")
                EXO_2_TP_side_map = [4 2 5 6 1 3];
                TP_side_id = EXO_2_TP_side_map(EXO_sideID);
            end
        end
    end
end

