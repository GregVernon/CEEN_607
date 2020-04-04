classdef feSolveParameters
    %feSolveParameters Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Incremental Nonlinear Iterations
        StepSize = 1;
        MaxSteps = 1;
        
        % Nonlinear Newton Iterations
        MaxNonlinearIterations = 10;
        RelNonlinearTolerance = 1e-12;
    end
    
    methods
        function obj = feSolveParameters()
            %feSolveParameters Construct an instance of this class
            %   Detailed explanation goes here
        end
    end
end

