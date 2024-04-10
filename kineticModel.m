classdef kineticModel
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Type 
        qmax {mustBeNumeric}
        Ks {mustBeNumeric}
        Ki {mustBeNumeric}
    end
    
    methods
        function obj = kineticModel(varargin)
            if nargin==0
                obj.Type='default';
            else
                obj.Type=varargin{1};
            end
            switch obj.Type
               case 'monod'
                  obj.qmax = varargin{2};
                  obj.Ks = varargin{3};
               case 'haldane'
                  obj.qmax = varargin{2};
                  obj.Ks = varargin{3};
                  obj.Ki = varargin{4};
                case 'inhibit'
                  obj.qmax = varargin{2};
                  obj.Ki = varargin{3};
                case 'proportional'
                   obj.qmax=1;
               otherwise
                  obj.qmax=1;
            end
        end
        
        function q = rate(obj,concentracion)
            switch obj.Type
               case 'monod'
                  q = obj.qmax*concentracion./(concentracion+obj.Ks);
               case 'haldane'
                  q = obj.qmax*concentracion./(concentracion+obj.Ks+concentracion.^2/obj.Ki);
                case 'inhibit'
                  q = obj.qmax*obj.Ki./(concentracion+obj.Ki);
                case 'proportional'
                  q= obj.qmax*concentracion;
               otherwise
                  q=1;
            end        
        end
        
        function [qMax,concMax] = peak(obj)
            switch obj.Type
               case 'monod'
                  concMax = inf;
                  qMax = obj.qmax;
               case 'haldane'
                  concMax = (obj.Ks*obj.Ki)^0.5;
                  qMax = obj.rate(concMax);
                case 'inhibit'
                  concMax = 0;
                  qMax = obj.rate(concMax);
               otherwise
                  concMax = 0;
                  qMax = obj.qmax;
            end     
        end
        
    end
end

