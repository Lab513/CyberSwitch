function draw_phase_space(inputs_time,p, LacI_vector,TetR_vector,varargin)

    ip = inputParser();
    ip.addParameter('BKGD_magnification_factor',3,@isscalar);
    ip.addParameter('BKGD_image','intensity',@(x) any(strcmp(x,{'intensity','orientation','rfp-gfp','none'})))
    ip.addParameter('quiver','on',@(x) any(strcmp(x,{'on','off'})))
    ip.addParameter('lines','off',@(x) any(strcmp(x,{'on','off'})))
    ip.parse(varargin{:});
    if ~strcmp(ip.Results.BKGD_image,'none')
        BMF = ip.Results.BKGD_magnification_factor;
        quivercolor = [1 1 1];
    else
        BMF = 1;
        quivercolor = [0 0 0];
    end

    % If we are using a static input (i.e. no averaging of the phase space,
    % just phase space for static concentratiosn of atc and iptg):
    if isvector(inputs_time)
        if isrow(inputs_time)
            inputs_time = inputs_time';
        end
        inputs_time(3) = 1; % Any value will do
    end

    % "Magnification" for a more precise heatmap background
    magnify = @(x) interp1(x,1:1/BMF:length(x));
    LacI_vec_plus = magnify(LacI_vector);
    TetR_vec_plus = magnify(TetR_vector);

    % Initiate derivatives arrays:
    dLacI=zeros(length(LacI_vec_plus),length(TetR_vec_plus));
    dTetR=zeros(length(LacI_vec_plus),length(TetR_vec_plus));
    
    for tim=1:size(inputs_time,2)
        % Get phase space for these copncentrations:
        [dydt] = get_phase_space(inputs_time(1:2,tim), p, LacI_vec_plus, TetR_vec_plus);
        % Time averaging: (nothing changes in static concentrations)
        dLacI = dLacI + dydt(:,:,3)*((inputs_time(3,tim)./sum(inputs_time(3,:))));
        dTetR = dTetR + dydt(:,:,4)*((inputs_time(3,tim)./sum(inputs_time(3,:))));
    end

    % Plot background image
    plot_bkgd_image(LacI_vector,TetR_vector,dLacI,dTetR,ip.Results.BKGD_image);
    hold on

    
    % Plot quiver
    if strcmp(ip.Results.quiver,'on')
        quiver(     repmat(LacI_vector',size(TetR_vector)), ...
                    repmat(TetR_vector,size(LacI_vector')), ...
                    dLacI(1:BMF:end,1:BMF:end),             ...
                    dTetR(1:BMF:end,1:BMF:end),             ...
                    2,'Color',quivercolor,'LineWidth',2); %the number 3 is just to scale the arrows...
    end

    % Cosmetics
    set(gca,'YDir','normal')
    xlabel('LacI (copy numbers)')
    ylabel('TetR (copy numbers)')
    ylim([min(TetR_vector) max(TetR_vector)])
    xlim([min(LacI_vector) max(LacI_vector)])
    set(gca,'TickDir','out')
    set(gca,'box','off')
end

function [dydt] = get_phase_space(inputs, p, LacI_vector, TetR_vector)

    % Compute TetR QSSA mrna levels
    mRNAt = (p.crt ...
            + p.cit * hill_func(    LacI_vector * hill_func(    inputs(2), ...
                                                                p.kiptg, ...
                                                                p.miptg ...
                                                            ), ...
                                    p.k_l, ...
                                    p.nl ...
                                ) ...
            )...
            ./p.delta_mrnat;

    % Compute LacI QSSA mrna levels
    mRNAl = (p.crl ...
            + p.cil * hill_func(    TetR_vector * hill_func(    inputs(1), ...
                                                                p.katc, ...
                                                                p.matc ...
                                                            ), ...
                                    p.k_t, ...
                                    p.nt ...
                                )...
            ) ...
            ./p.delta_mrnal;

    grid = cat(3,repmat(mRNAl,size(TetR_vector')), repmat(mRNAt',size(LacI_vector)), repmat(LacI_vector',size(TetR_vector)), repmat(TetR_vector,size(LacI_vector')));
    dydt = toggle_derivative_sim(1,grid,inputs,p);
end

function plot_bkgd_image(LacI_vector,TetR_vector,dLacI,dTetR,BKGD_image)

    switch BKGD_image
        case 'intensity'
            BI = log10(sqrt(dLacI.^2 + dTetR.^2))';
            imagesc(LacI_vector([1,end]),TetR_vector([1,end]),BI);
            colormap('jet');
            c = colorbar('Ticks',-10:10,'TickLabels',10.^(-10:10));
            c.Label.String = 'Vector field intensity (prot/min)';
        case 'orientation' % TODO
            dTetRt = dTetR';
            dLacIt = dLacI';
            ang = atan2(dTetRt(:),dLacIt(:));
            BI = reshape(ang,length(dLacI),length(dTetR));
            imagesc(LacI_vector([1,end]),TetR_vector([1,end]),BI);
            colormap('hsv');
            c = colorbar('Ticks',-pi:(pi/4):pi,'TickLabels',{'-\pi' '-3\pi/4' '-\pi/2' '-\pi/4' '0' '\pi/4' '\pi/2' '3\pi/4' '\pi' });
            c.Label.String = 'Vector field orientation';
        case 'rfp-gfp'
            BI(:,:,1) = repmat(LacI_vector',size(TetR_vector))./max(LacI_vector);
            BI(:,:,2) = repmat(TetR_vector,size(LacI_vector'))./max(TetR_vector);
            BI(:,:,3) = 0;
            imagesc(LacI_vector([1,end]),TetR_vector([1,end]),BI);
        case 'none'
            return;
    end
end