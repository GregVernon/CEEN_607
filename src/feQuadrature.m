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
    
    methods (Static)
        function [P,W] = assembleQuadratureRule(nDim, NumPoints, qID)
            if nDim == 1
                P = cell(NumPoints,1);
                W = cell(NumPoints,1);
                [P1,W1] = feQuadrature.GaussLegendre(NumPoints(1));
                qp = 0;
                for ii = 1:NumPoints(1)
                    qp = qp+1;
                    P{qp} = P1(ii);
                    W{qp} = W1(ii);
                end
            elseif nDim == 2
                P = cell(prod(NumPoints),1);
                W = cell(prod(NumPoints),1);
                if qID == 0
                    [P1,W1] = feQuadrature.GaussLegendre(NumPoints(1));
                    [P2,W2] = feQuadrature.GaussLegendre(NumPoints(2));
                    qp = 0;
                    for jj = 1:NumPoints(2)
                        for ii = 1:NumPoints(1)
                            qp = qp+1;
                            P{qp} = [P1(ii); P2(jj)];
                            W{qp} = W1(ii) * W2(jj);
                        end
                    end
                elseif qID == 1 % XMin side
                    [P2,W2] = feQuadrature.GaussLegendre(NumPoints);
                    qp = 0;
                    for ii = NumPoints:-1:1
                        qp = qp+1;
                        P{qp} = [-1; P2(ii)];
                        W{qp} = W2(ii);
                    end
                elseif qID == 2
                    [P2,W2] = feQuadrature.GaussLegendre(NumPoints);
                    qp = 0;
                    for ii = 1:NumPoints
                        qp = qp+1;
                        P{qp} = [1; P2(ii)];
                        W{qp} = W2(ii);
                    end
                elseif qID == 3 
                    [P1,W1] = feQuadrature.GaussLegendre(NumPoints);
                    qp = 0;
                    for ii = NumPoints:-1:1
                        qp = qp+1;
                        P{qp} = [P1(ii); -1];
                        W{qp} = W1(ii);
                    end
                elseif qID == 4
                    [P1,W1] = feQuadrature.GaussLegendre(NumPoints);
                    qp = 0;
                    for ii = 1:NumPoints
                        qp = qp+1;
                        P{qp} = [P1(ii); 1];
                        W{qp} = W1(ii);
                    end
                end
            elseif nDim == 3
                P = cell(prod(NumPoints),1);
                W = cell(prod(NumPoints),1);
                [P1,W1] = feQuadrature.GaussLegendre(NumPoints(1));
                [P2,W2] = feQuadrature.GaussLegendre(NumPoints(2));
                [P3,W3] = feQuadrature.GaussLegendre(NumPoints(3));
                qp = 0;
                for kk = 1:NumPoints(3)
                    for jj = 1:NumPoints(2)
                        for ii = 1:NumPoints(1)
                            qp = qp+1;
                            P{qp} = [P1(ii) P2(jj) P3(kk)];
                            W{qp} = W1(ii) * W2(jj) * W3(kk);
                        end
                    end
                end
            end
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
                    P = [-1/sqrt(3); +1/sqrt(3)];
                    W = [1; 1];
            end
        end
    end
end

