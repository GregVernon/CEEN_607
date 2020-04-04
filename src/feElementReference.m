classdef feElementReference
    %feElementConfiguration Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Descriptors
        Degree
        Dimension
        GlobalID
        Type
                
        % Children
        Nodes
        Quadrature
        
        % Physics
        MaterialConstitutiveMatrix
        
        % Boundary Conditions
        DirichletConditions
        
        % Load Conditions
        BodyForce
        NodeForce
        SurfacePressure
        
        % Linear Algebra Stuff
        StiffnessMatrix
        InternalForceVector
        ExternalForceVector
        
        % Solution Container
        Solution
    end
    
    methods
        % Constructors
        function obj = feElementReference()
            %feElementConfiguration Construct an instance of this class
        end
    end
end

