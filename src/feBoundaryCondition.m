classdef feBoundaryCondition
    %feBoundaryCondition Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Type
        ElementSet feElementSet
        SurfaceSet feSurfaceSet
        NodeSet    feNodeSet
        DOF
        Value
    end
    
    methods
        function obj = feBoundaryCondition()
            %feBoundaryCondition Construct an instance of this class
            %   Detailed explanation goes here
        end
    end
end

