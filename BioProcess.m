classdef BioProcess < dynamicprops 
    % Bioprocess model class
    %   Class to define bioprocess dynamic model properties and perform
    %   calculations.
    %   Instances must be created specifying number of states and
    %   reactions.
    
    properties
        NumberStates {mustBeNumeric}
        NumberReactions {mustBeNumeric}
        YieldMatrix {mustBeNumeric}
        RateArray cell
        Names cell
    end
    
    methods
        function obj=BioProcess(numStates,numReactions)  
            % Constructor
            obj.NumberStates=numStates;
            obj.NumberReactions=numReactions;
            obj.YieldMatrix=nan(numStates,numReactions);
            C = cell(numReactions,numStates);
            for idx=1:numReactions
                for jdx=1:numStates
                    C{idx,jdx}=kineticModel();
                end
            end
            obj.RateArray = C;
%             obj.Names= cell(1,numStates);
        end
        
        function set.YieldMatrix(obj,value)
            %set YieldMatrix Property. 
            % Matrix must be NumberStates x NumberReactions
            if size(value)== [obj.NumberStates,obj.NumberReactions]
                obj.YieldMatrix=value;           
            else
                error('Incorrect dimensions, matrix must be %d by %d',obj.NumberStates,obj.NumberReactions);
            end
        end
        
        function set.RateArray(obj,value)
            % set the kinetic rate factors
            % the cell array RateArray contains the kinetic models factors 
            % for each reaction rate with respect to each state variable.
            % Each element of the aaray must belong to the kineticModel
            % class
            if iscell(value)
                dim=size(value);
                errores=0;
                for idx=1:dim(1)*dim(2)
                    if ~isequal(class(value{idx}),'kineticModel')                     
                       errores=errores+1;
                    end
                end
                if errores~=0
                    error('One or more inserted elements do not belong to the class kineticModel'); 
                else
                    obj.RateArray=value; 
                end
            else
                if isequal(class(value),'kineticModel')
                    obj.RateArray=value; 
                else
                    error('The inserted element does not belong to the class kineticModel');                  
                end
            end
        end     
        
        function set.Names(obj,value)
            % set the names of the state variables
            % the names are stored in a cell array
            % must be strings
            if iscell(value)
                if size(value)==[1,obj.NumberStates]
                    errores=0;
                    for idx=1:obj.NumberStates
                        if ~isstring(value{idx}) && ~ischar(value{idx})                       
                           errores=errores+1;
                        end
                    end
                    if errores~=0
                        error('One or more inserted elements is not a string'); 
                    else
                        obj.Names=value; 
                    end
                else
                    error('The size of the cells do not match. Must have %d names',obj.NumberStates);
                end
            else
                if ~isstring(value{idx}) || ~ischar(value{idx})
                    error('The inserted elementis not a char or string'); 
                else
                    obj.Names=value; 
                end
            end
            
        end
       
        function dxidt = stateSpace(obj,t,xi,xi_in,D)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            dxidt = obj.YieldMatrix*obj.rates(xi)+D*(xi_in-xi);
        end
        
        function r = rates(obj,xi)
            r=ones(obj.NumberReactions,1);
            for idx=1:obj.NumberReactions
                for jdx=1:obj.NumberStates
                    modelo = obj.RateArray{idx,jdx};
                    r(idx)= r(idx).*modelo.rate(xi(jdx));
                end
            end
        end
        
        function [t,xi,varargout] = euler(obj,dt,tspan,xi0,xi_in,Dfeed,varargin)
            args=length(varargin);
            t=tspan(1):dt:tspan(2);
            N=length(t);
            xi=nan(obj.NumberStates,N);
            xi(:,1)=xi0;
            if args>=2
                V=nan(1,N);
                V(1)=varargin{2};
            end
            
            for idx=2:N
                if args>0
                    tFeed = varargin{1};
                    if idx*dt>=tFeed
                        D=Dfeed;
                    else
                        D=0;
                    end
                end
                xi(:,idx)=xi(:,idx-1)+dt*obj.stateSpace(t,xi(:,idx-1),xi_in,D);
                if args>=2
                    V(idx)=V(idx-1)+dt*D*V(idx-1);
                end
            end
            
            t=t';
            xi=xi';
            if args>=2
                    varargout{1}=V';
            end
            
        end
        
    end
end

