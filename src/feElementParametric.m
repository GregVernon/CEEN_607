classdef feElementParametric
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
        
        % Linear Algebra Stuff
        StiffnessMatrix
        InternalForceVector
        ExternalForceVector
    end
    
    methods
        % Constructors
        function obj = feElementParametric()
            %feElementConfiguration Construct an instance of this class
        end
    end
    
    methods (Static)
        % Nodal Information
        function nodalCoords = computeNodeCoordinates(ElementDegree, loc_nod_id)
            if     length(ElementDegree) == 1
                x = linspace(-1,1,ElementDegree(1));
                nodalCoords = [x(loc_nod_id)];
            elseif length(ElementDegree) == 2
                nID = 0;
                x = linspace(-1,1,ElementDegree(1)+1);
                y = linspace(-1,1,ElementDegree(2)+1);
                for ny = 1:ElementDegree(2)+1
                    for nx = 1:ElementDegree(1)+1
                        nID = nID + 1;
                        if nID == loc_nod_id
                            nodalCoords = [x(nx); y(ny)];
                            return
                        end
                    end
                end
            elseif length(ElementDegree) == 3
                nID = 0;
                x = linspace(-1,1,ElementDegree(1));
                y = linspace(-1,1,ElementDegree(2));
                z = linspace(-1,1,ElementDegree(3));
                for nz = 1:ElementDegree(3)+1
                    for ny = 1:ElementDegree(2)+1
                        for nx = 1:ElementDegree(1)+1
                            nID = nID + 1;
                            if nID == loc_nod_id
                                nodalCoords = [x(nx); y(ny); z(nz)];
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end

