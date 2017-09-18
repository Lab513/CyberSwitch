% code for cusp catastrophe plots from JB
function [MSS, BSS, USS] = get_intersections_line(p)
%Greg: updated inputs order

    IPTG = linspace(0,.45,200);
    ATC = 45-(IPTG*100);

    MSS = zeros(1,4);
    BSS = zeros(1,4);
    USS = zeros(1,4);

    counter = 0;

    for inp1 = 1:length(IPTG)
        inp2 = inp1;
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
    %         xlabel('LacI (copy numbers)')
    %         ylabel('TetR (copy numbers)')

            %     plot vector field
            %     LacI_vector= linspace(0,2000,35);
            %     TetR_vector= linspace(0,2000,35);
            %     for i=1:length(LacI_vector)
            %         mRNAt(i) = (p.crt + p.cit./(1 + (LacI_vector(i)/p.k_l.*(1./(1 + inputs(1)./(p.kiptg*LacI_vector(i)))).^p.miptg).^p.nl))./p.delta_mrnal;
            %     end
            %     for j=1:length(TetR_vector)
            %         mRNAl(j) = (p.crl + p.cil./(1 + (TetR_vector(j)/p.k_t.*(1./(1 + inputs(2)./(p.katc*TetR_vector(j)))).^p.matc).^p.nt))./p.delta_mrnat;
            %     end
            %     for i= 1:length(LacI_vector)
            %         for j= 1:length(TetR_vector)
            %             grid(i,j,:)= [LacI_vector(i) TetR_vector(j)];
            %             dydt= toggle_derivative2(0,[mRNAl(j) mRNAt(i) LacI_vector(i) TetR_vector(j)],inputs,p,set);
            %             dLacI(i,j)= dydt(3);
            %             dTetR(i,j)= dydt(4);
            %         end
            %     end
            %     quiver(grid(:,:,1),grid(:,:,2),dLacI,dTetR,3); %the number 3 is just to scale the arrows...

        fprintf('%.02f %% Done \n', 100*counter/(length(IPTG)*length(ATC)));
    end
end