function plot_mean_optimised_results( p , data , varargin)
    
ip = inputParser;
addParameter(ip,'RunSSA',false,@islogical);
ip.parse(varargin{:});

    for ds=1:numel(data)

        timegfp = data(ds).val.timegfp;
        timerfp = data(ds).val.timerfp;
        gfpMothers = data(ds).val.gfpMothers;
        rfpMothers = data(ds).val.rfpMothers;
        inputs =  data(ds).ind.inputs;
        inputsbefore = data(ds).ind.inputsbefore;
        
        fprintf('Simulating dataset #%d\n',ds) 
        disp('here: used observed initial conditions')
        Y0= startingpoint_value(p, inputsbefore, data(ds).val);
        %Y0=startingpoint_ode( p ,inputsbefore, data(ds).val);
        
        [~,~,LacI_opt,TetR_opt,iptg,atc] = generate_data(p,Y0,timerfp/60,inputs);
        if ip.Results.RunSSA
            parfor indCell = 1:12
                [~,~,LacIs(indCell,:),TetRs(indCell,:),~,~] = generate_data(p,Y0,timerfp/60,inputs,'simulation_method','SSA');
            end
        else
            LacIs = [];
            TetRs = [];
        end

        % Plot original data:
        subplot(numel(data),3,ds*3-1)
        hold on
        plot(timerfp/3600,mean(rfpMothers,2),'r','LineWidth',2);
        plot(timerfp/3600,rfpMothers,'r');
        plot(timerfp/3600,mean(gfpMothers,2),'g','LineWidth',2);
        plot(timerfp/3600,gfpMothers,'g');
        axis([0 length(inputs)/60 0 4.5e3]);
        
        % Plot optimized model
        subplot(numel(data),3,ds*3)
        hold on
        plot(timerfp/3600,TetR_opt,'g','linewidth',2);
        plot(timerfp/3600,TetRs,'Color','g');
        plot(timerfp/3600,LacI_opt,'r','linewidth',2);
        plot(timerfp/3600,LacIs,'r');
        
        clearvars LacIs TetRs

        axis([0 length(inputs)/60 0 4.5e3]);
%         if isfield(data(ds),'name')
%             title(data(ds).name,'Interpreter','None')
%         end
     
        subplot(numel(data),3,ds*3-2);
        [AX,H1,H2] = plotyy((1:length(inputs))/60,inputs(2,:)',(1:length(inputs))/60,inputs(1,:)',@stairs,@stairs);
        set(H1,'Color','c')
        set(H1,'Linewidth',2)
        set(H2,'Color','m')
        set(H2,'Linewidth',2)
        set(AX(1),'Ycolor','c')
        set(get(AX(1),'Ylabel'),'String','IPTG') 
        set(AX(1),'Ylim',[-0.05 1.05])
        set(AX(2),'Ycolor','m')
        set(get(AX(2),'Ylabel'),'String','aTC') 
        set(AX(2),'Ylim',[-5 105])
        xlabel('time (h)');
        set(AX(1),'Xlim',[0 length(inputs)/60])
        set(AX(2),'Xlim',[0 length(inputs)/60])
        hold on; plot(timerfp/3600,iptg,'c--','Linewidth',2);
        hold on; plot(timerfp/3600,atc/100,'m--','Linewidth',2);
        drawnow
    end
end