% code for cusp catastrophe plots from JB
function [MSS, BSS, USS] = get_intersections_jb(p)
%Greg: updated inputs order

    IPTG = linspace(0,1,100);
    ATC = linspace(0,100,100);

    MSS = zeros(1,4);
    BSS = zeros(1,4);
    USS = zeros(1,4);

    counter = 0;

    for inp1 = 1:length(IPTG)
        for inp2 = 1:length(ATC)
            counter = counter + 1;
            inputs(2) = IPTG(inp1);
            inputs(1) = ATC(inp2);
            LacI_vector= linspace(0,4000,300);
            TetR_vector= linspace(0,1500,300);
            [X0, Y0] = get_intersections(inputs,p, LacI_vector,TetR_vector);
            if length(X0)==1
    %             plot(X0,Y0,'ko'); hold on;
                MSS(end+1,:) = [ATC(inp2) IPTG(inp1) X0(1) Y0(1)];
            else
    %             plot(X0(1),Y0(1),'mo'); hold on;
                BSS(end+1,:) = [ATC(inp2) IPTG(inp1) X0(1) Y0(1)];
    %             plot(X0(2),Y0(2),'ro'); hold on;
                USS(end+1,:) = [ATC(inp2) IPTG(inp1) X0(2) Y0(2)];
    %             plot(X0(3),Y0(3),'mo'); hold on;
                BSS(end+1,:) = [ATC(inp2) IPTG(inp1) X0(3) Y0(3)];
            end

        end
        fprintf('%.02f %% Done \n', 100*counter/(length(IPTG)*length(ATC)));
    end
end