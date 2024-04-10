classdef BioProcess < dynamicprops 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NumberStates {mustBeNumeric}
        NumberReactions {mustBeNumeric}
        YieldMatrix {mustBeNumeric}
        RateArray cell
    end
    
    methods
        function obj=BioProcess(numStates,numReactions)  
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
        end
        
        function obj = set.YieldMatrix(obj,value)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.YieldMatrix=value;           
        end
        
        function obj = set.RateArray(obj,value)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here 
            obj.RateArray=value; 
        end
        
        function obj = get(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function r = reactionRates(obj,xi)
            r=nan(obj.NumberReactions,1);
            
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
                    r(idx)= r(idx)*modelo.rate(xi(jdx));
                end
            end
        end
        
    end
end

