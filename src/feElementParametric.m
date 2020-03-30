classdef feElementParametric
    %feElementConfiguration Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Descriptors
        ElementType
        ElementDegree
        ElementGlobalID
        
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
        function nodalCoords = getNodeCoordinates(ElementDegree, loc_nod_id)
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
    
    methods
        % Construct Children
        function obj = createChildrenNodes(obj)
            nDim = length(obj.ElementDegree);
            nNodes = obj.ElementDegree+1;
            NumTotalNodes = prod(nNodes);
            obj.Nodes = repmat(feNode("Parametric"),NumTotalNodes,1);
            for n = 1:NumTotalNodes
                obj.Nodes(n).ConfigurationType = "Parametric";
                obj.Nodes(n).Coordinates = feElementParametric.getNodeCoordinates(obj.ElementDegree, n);
                obj.Nodes(n).Dimension = nDim;
            end
        end
        
        function obj = createQuadrature(obj)
            nDim = length(obj.ElementDegree);
            nNodes = obj.ElementDegree+1;
            NumTotalNodes = prod(nNodes);
            if     nDim == 1
                % Single quadrature entity on a bar element domain
                obj.Quadrature = feQuadrature("Parametric");
            elseif nDim == 2
                % Quadrature entities on body + 4 sides
                obj.Quadrature = repmat(feQuadrature("Parametric"),5,1);
            elseif nDim == 3
                % Quadrature entities on body + 6 sides
                obj.Quadrature = repmat(feQuadrature("Parametric"),7,1);
            end
        end
    end
end

