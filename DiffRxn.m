classdef DiffRxn
    %DIFFRXN - Structure to test VPAM with Diffusion Reaction Problem
    %   Holds the data associated with the diffusion-reaction problem
    %   and contains methods for estimating non-linear residuals
    %   and other information associated.
    
    properties
        %   Data objects
        DiffusionCoef
        ReactionCoef
        BoundaryVal
        Length
        BasisSize
    end
    
    methods
        % Function objects
        function obj = DiffRxn()
            %DIFFRXN Construct an instance of this class
            %   Set default values for all properties 
            obj.DiffusionCoef = 1;
            obj.ReactionCoef = 1;
            obj.BoundaryVal = 1;
            obj.Length = 1;
            obj.BasisSize = 3;
        end
        
        function F = Residual(obj, x)
            %Residual - Function to compute VPAM residuals from data
            %   Takes data in object and computes all integrals for
            %   all functions in residual.
            
            % Loop to compute boundary conditions 
            res0 = 0; res1 = 0;
            for i=1:obj.BasisSize
                res0 = res0 + (x(i+2)*PolyBasis1D(i, 0.0));
                res1 = res1 + (x(i+2)*FirstDerivativePolyBasis1D(i, obj.Length));
            end
            F(1,1) = res0 - obj.BoundaryVal;
            F(2,1) = res1;
            
            % Loop for all other residuals 
            for i=1:obj.BasisSize
                Lap_sum = 0; Over_sum = 0;
                for j=1:obj.BasisSize
                    Lap_sum = Lap_sum + ( x(j+2)*LapIntegral(i,j,0,obj.Length) );
                    Over_sum = Over_sum + ( x(j+2)*OverlapIntegral(i,j,0,obj.Length) );
                end
                F(i+2,1) = Lap_sum - (obj.ReactionCoef*Over_sum/obj.DiffusionCoef) + (x(1)*PolyBasis1D(i,0)) + (x(2)*FirstDerivativePolyBasis1D(i,obj.Length));
            end
        end
        
        function [u, z] = Evaluate(obj, x, N)
        %Evaluate - Evaluate the solution for u given vector x
        %   Note: x(1) and x(2) are the lambda's (don't use)
        %   Returns a vector with solutions at different points
            switch nargin
                case 3
                    % nothing
                case 2
                    N = 20;
                case 1
                    N = 20;
                    x = zeros(obj.BasisSize+2,1);
                otherwise
                    N = 20;
                    x = zeros(obj.BasisSize+2,1);
            end 
            
            c = zeros(obj.BasisSize,1);
            for i=1:obj.BasisSize
                c(i) = x(i+2);
            end    
            u = zeros(N+1,1);
            z = zeros(N+1,1);
            pos = 0;
            for i=1:N+1
                z(i) = pos;
                u(i) = EvalPoly(c,pos);
                pos = pos + (obj.Length/N);
            end
            
        end
        
        function [u, z] = ExactSoln(obj, N)
        %Evaluate - Evaluate the solution for u given vector x
        %   Note: x(1) and x(2) are the lambda's (don't use)
        %   Returns a vector with solutions at different points
            switch nargin
                case 2
                    % nothing
                case 1
                    N = 20;
                otherwise
                    N = 20;
            end 
                
            u = zeros(N+1,1);
            z = zeros(N+1,1);
            pos = 0;
            L = obj.Length;         
            k = obj.ReactionCoef;   
            D = obj.DiffusionCoef;  
            uo = obj.BoundaryVal;   
            for i=1:N+1
                z(i) = pos;
                u(i) = (uo/(1+exp(2*sqrt(k/D)*L)))*(exp(sqrt(k/D)*((2*L)-pos))+exp(sqrt(k/D)*pos));
                pos = pos + (obj.Length/N);
            end
            
        end
        
    end

end

