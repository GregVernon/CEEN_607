classdef feMesh
    %feMesh Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
    end
end

