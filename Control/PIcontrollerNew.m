classdef PIcontrollerNew < handle
    properties
        P = 1/50; % Those are roughly estimated for both branch, but we should use finer parameters when instantiating
        I = 1/500;
        idelay = 0; % In minutes
        objective = [2000,Inf]; % Some random objective, should be changed depending on the branch
        capdecision = [0 50];
        levelzero = 10; % The level that is considered 'no induction' for the PWM (to be changed depending on the branch!)
    end
    
    methods
        function obj = PIcontrollerNew(varargin)
            if nargin >= 1
                obj.P = varargin{1};
            end
            if nargin >= 2
                obj.I = varargin{2};
            end
            if nargin >= 3
                obj.objective = varargin{3};
            end
            if nargin >= 4
                obj.capdecision = varargin{4};
            end
            if nargin >= 5
                obj.levelzero = varargin{5};
            end
            if nargin >= 6
                obj.idelay = varargin{6};
            end
        end
        
        function decision = decide(obj,values,tpoints)
            c = obj.capdecision;
            err = values - obj.reconstructobjective(tpoints)';
            
            
            try % The delay version
                start_int= find(tpoints>obj.idelay*60, 1, 'first');
                if ~isempty(start_int) && numel(tpoints) - start_int > 1
                    interr = trapz(tpoints(start_int:end)./60,err(start_int:end));
                else
                    interr = 0;
                end
            catch % If doesn't work, fall back to normal version
                disp('catch here');
                if numel(tpoints) > 1
                    interr = trapz(tpoints./60,err);
                else
                    interr = 0;
                end
            end
            PropD = obj.P*err(end);
            IntegD = obj.I*interr;
            decision = obj.levelzero -(PropD + IntegD); % I use timepoints in minutes here! Remember when setting the I parameter
            decisionO = round(decision);
            decision(decision<c(1)) = c(1);
            decision(decision>c(2)) = c(2);
            decision = round(decision);
            try
                %disp(['[PIcontroller] ' datestr(now) ' - error: ' num2str(round(err(end))) ', int. error: ' num2str(round(interr)) ' (' num2str(numel(tpoints)) ' execs), decision: ' num2str(decisionO) ' (P=' num2str(round(PropD)) ', I=' num2str(round(IntegD)) ') --> ' num2str(decision)])
            catch 
            end
        end
        
        function objrec = reconstructobjective(obj,tpoints)
            % objrec = interp1([-1e6 -1e-6 0 obj.objective(:,2)'],[] %WTF?
            % 'previous' and 'next' give me terribly wrong results!
            % I have to do this for loop instead:
            
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