classdef feNode
    %feNode Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Topology
        ChildDOF_ID = missing
        
        % Descriptors
        ConfigurationType = missing
        Coordinates = missing
        Dimension = missing
        GlobalID = missing
    end
    
    methods
        % Constructors
        function obj = feNode(ConfigurationType)
            %feNode Construct an instance of this class
            %   Detailed explanation goes here
            obj.ConfigurationType = ConfigurationType;
        end
    end
    
    methods

    end
end

