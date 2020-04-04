classdef feElementSet
    %feElementSet Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        ElementID
        ElementType
        NodeConnectivity
    end
    
    methods
        % Constructors
        function obj = feElementSet()
            %feElementSet Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = build_from_exodus(obj, EXO, setName)
            % Exodus.element_blocks{:,set_id}
                % 1 - name
                % 2 - id
                % 3 - topology type
                % 4 - connectivity
            set_id = find(strcmpi(setName, EXO.element_blocks(1,:)));
            obj.Name = EXO.element_blocks{1,set_id};
            obj.ElementID = EXO.element_blocks{2,set_id};
            obj.ElementType = EXO.element_blocks{3,set_id};
            obj.NodeConnectivity = EXO.element_blocks{4,set_id};
        end
    end
end

