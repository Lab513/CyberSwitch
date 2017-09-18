classdef BBcontrollerNew < handle
    properties
        objective = [2000,Inf]; % Some random objective, should be changed depending on the branch
        capdecision = [0 50];
    end
    
    methods
        function obj = BBcontrollerNew(varargin)
            if nargin >= 1
                obj.objective = varargin{1};
            end
            if nargin >= 2
                obj.capdecision = varargin{2};
            end
        end
        
        function decision = decide(obj,values,tpoints)
            c = obj.capdecision;
            err = values - obj.reconstructobjective(tpoints)';
            
            if err(end) > 0
                decision = c(1);
            else
                decision = c(2);
            end
            
            try
                %disp(['[PIcontroller] ' datestr(now) ' - error: ' num2str(round(err(end))) ', int. error: ' num2str(round(interr)) ' (' num2str(numel(tpoints)) ' execs), decision: ' num2str(decisionO) ' (P=' num2str(round(PropD)) ', I=' num2str(round(IntegD)) ') --> ' num2str(decision)])
            catch 
            end
        end
        
        function objrec = reconstructobjective(obj,tpoints)
            
            lastobj = tpoints(1);
            for ind1 = 1:size(obj.objective)
                objrec(tpoints < obj.objective(ind1,2) & tpoints >= lastobj) = obj.objective(ind1,1);
                lastobj = obj.objective(ind1,2);
            end
            objrec(tpoints >= lastobj) = obj.objective(ind1,1); % If there are points beyond the last objective given
        end  
        
        
        function visp(obj,msg,varargin)

            % If no verbosity specified, verbosity is one
            if nargin == 2
                verbosity = 1;
            else
                verbosity = varargin{1};
            end
            
            % Get level of verbosity of object:
            verb = obj.verbosity;
            
            % Cat message and output:
            msg = ['[' class(obj) '] ' datestr(now) ' - '   msg '\n'];
            if verbosity < 0
                error(msg);
            end
            if verb >= verbosity
                fprintf(msg);
            end
           
        end
        
    end
    
end