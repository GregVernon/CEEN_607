classdef feQuadrature
    %feQuadrature Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Descriptors
        IntegralType
        Scheme
        % IntegralType: ALL
        NumPoints
        Coordinates
        Weights
        BasisFunction
        IntegralScaling
        % IntegralType: SOLID
        JacobianMatrix
        StrainDisplacementMatrix
        % IntegralType: BOUNDARY
        NormalVector
        TangentVector
    end
    
    methods
        % Constructor        
        function obj = feQuadrature(IntegralType, Scheme)
            %feQuadrature Construct an instance of this class
                % Builds the class with a quadrature scheme with nPts
                % If initZeros is true, just preallocate arrays (zeros)
            obj.IntegralType = IntegralType;
            obj.Scheme = Scheme;
        end
    end
    
    methods
        function obj = buildParametricCoordinates(obj, Scheme, nPts)
            if nargin > 1
                obj.Scheme = Scheme;
                obj.NumPoints = nPts;
            end
            nDim = length(obj.NumPoints);
            if nDim == 1
                P = cell(obj.NumPoints);
                W = cell(obj.NumPoints);
                [P1,W1] = obj.GaussLegendre(obj.NumPoints(1));
                qp = 0;
                for ii = 1:obj.NumPoints(1)
                    qp = qp+1;
                    P{qp} = P1(ii);
                    W{qp} = W1(ii);
                end
            elseif nDim == 2
                P = cell(prod(obj.NumPoints));
                W = cell(prod(obj.NumPoints));
                [P1,W1] = obj.GaussLegendre(obj.NumPoints(1));
                [P2,W2] = obj.GaussLegendre(obj.NumPoints(2));
                qp = 0;
                for jj = 1:obj.NumPoints(2)
                    for ii = 1:obj.NumPoints(1)
                        qp = qp+1;
                        P{qp} = [P1(ii) P2(jj)];
                        W{qp} = W1(ii) * W2(jj);
                    end
                end
            elseif nDim == 3
                P = cell(prod(obj.NumPoints));
                W = cell(prod(obj.NumPoints));
                [P1,W1] = obj.GaussLegendre(obj.NumPoints(1));
                [P2,W2] = obj.GaussLegendre(obj.NumPoints(2));
                [P3,W3] = obj.GaussLegendre(obj.NumPoints(3));
                qp = 0;
                for kk = 1:obj.NumPoints(3)
                    for jj = 1:obj.NumPoints(2)
                        for ii = 1:obj.NumPoints(1)
                            qp = qp+1;
                            P{qp} = [P1(ii) P2(jj) P3(kk)];
                            W{qp} = W1(ii) * W2(jj) * W3(kk);
                        end
                    end
                end
            end
            obj.Coordinates = P;
            obj.Weights = W;
        end
    end
    
    methods (Static)
        % Quadrature schemes
        function [P,W] = GaussLegendre(nPts)
            %GaussLegendre Return nPts Rule for a 1D Gauss-Legendre Scheme
            %   Detailed explanation goes here
            switch nPts
                case 1
                    P = 0;
                    W = 2;
                case 2
                    P = [-1/sqrt(3) +1/sqrt(3)];
                    W = [1 1];
            end
        end
    end
end

