function [Y,xY] = my_spm_regions(SPM, my_xyz, my_name, my_con_nr, my_radius)

%% my_xyz == xyz of VOI centre
%% my_name == name of VOI
%% my_con_nr == contrast used for adjustment an 1 needs to be added
%% my_radius == radius of VOI when sphere is used

% VOI time-series extraction of adjusted data (& local eigenimage analysis)
% FORMAT [Y xY] = spm_regions(xSPM,SPM,hReg,[xY]);
%
% xSPM   - structure containing specific SPM, distribution & filtering details
% SPM    - structure containing generic analysis details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Y      - first scaled eigenvariate of VOI {i.e. weighted mean}
% xY     - VOI structure
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.Ic           - contrast used to adjust data (0 - no adjustment)
%       xY.Sess         - session index
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.XYZmm        - Co-ordinates of VOI voxels {mm}
%       xY.y            - [whitened and filtered] voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues
%       xY.X0           - [whitened] confounds (including drift terms)
%
% Y and xY are also saved in VOI_*.mat in the SPM working directory
%
% (See spm_getSPM for details on the SPM,VOL, xX & xSDM structures.)
%
%_______________________________________________________________________
%
% spm_regions extracts a representative time course from voxel data
% in terms of the first eigenvariate of the filtered and adjusted
% response in all suprathreshold voxels within a specified VOI
% centered on the current MIP cursor location.
%
% If temporal filtering has been specified, then the data will be
% filtered.  Similarly for whitening. Adjustment is with respect to
% the null space of a selected contrast, or can be omitted.
%
% For a VOI of radius 0, the [adjusted] voxel time-series is
% returned, and scaled to have a 2-norm or 1. The actual [adjusted]
% voxel time series can be extracted from xY.y, and will be
% the same as the [adjusted] data returned by the plotting routine
% (spm_graph.m) for the same contrast.
%_______________________________________________________________________
% @(#)spm_regions.m	2.18 Karl Friston 03/11/03

%%% Get all voxels suprathreshold
xSPM.XYZ=SPM.xVol.XYZ;
XYZtemp = [SPM.xVol.XYZ; ones(1,length(SPM.xVol.XYZ))];
XYZmm =SPM.xVol.M * XYZtemp;
xSPM.XYZmm =XYZmm(1:3,:);


%%%% Konstanten fuer VOI
myVOI_spec = 'sphere'


% get figure handles
%-----------------------------------------------------------------------
% Finter = spm_figure('GetWin','Interactive');
% Fgraph = spm_figure('GetWin','Graphics');
% header = get(Finter,'Name');
% set(Finter,'Name','VOI time-series extraction')
% try
% 	xY;
% catch
xY = {};
% end

% %-Find nearest voxel [Euclidean distance] in point list in Y.mad
% %-----------------------------------------------------------------------
% if ~length(xSPM.XYZmm)
% 	spm('alert!','No suprathreshold voxels!',mfilename,0);
% 	Y = []; xY = [];
% 	return
% end
% try
% 	xyz     = xY.xyz;
% catch
% 	[xyz,i] = spm_XYZreg('NearestXYZ',...
% 		      spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
% 	xY.xyz  = xyz;
% end
%
% % and update GUI location
% %-----------------------------------------------------------------------
% spm_XYZreg('SetCoords',xyz,hReg);


xyz=my_xyz
xY.xyz  = xyz;
%

%-Get adjustment options and VOI name
%-----------------------------------------------------------------------
% spm_input(sprintf('at [%3.0f %3.0f %3.0f]',xY.xyz),1,'d',...
% 	'VOI time-series extraction')

if ~isfield(xY,'name')
    xY.name = my_name;
end

if ~isfield(xY,'Ic')
    q     = 0;
    Con   = {'<don''t adjust>'};
    for i = 1:length(SPM.xCon)
        if strcmp(SPM.xCon(i).STAT,'F')
            q(end + 1) = i;
            Con{end + 1} = SPM.xCon(i).name;
        end
    end
    i     = my_con_nr;   %%%%number of contrast
    xY.Ic = q(i);
end

%-if fMRI data get sessions and filtering options
%-----------------------------------------------------------------------
for s = 1:length(SPM.Sess)
    xY.Sess   = s;
    Sessname = {['sess' num2str(s)]}


    %-Specify VOI
    %-----------------------------------------------------------------------
    if ~isfield(xY,'def')
        xY.def    = myVOI_spec;
    end
    Q       = ones(1,size(xSPM.XYZmm,2));


    switch xY.def

        case 'sphere'
            %---------------------------------------------------------------
            if ~isfield(xY,'spec')
                xY.spec = my_radius;
            end
            d     = [xSPM.XYZmm(1,:) - xyz(1);
                xSPM.XYZmm(2,:) - xyz(2);
                xSPM.XYZmm(3,:) - xyz(3)];
            Q     = find(sum(d.^2) <= xY.spec^2);

        case 'box'
            %---------------------------------------------------------------
            if ~isfield(xY,'spec')
                xY.spec = spm_input('box dimensions [x y z] {mm}',...
                    '!+0','r','0 0 0',3);
            end
            Q     = find(all(abs(xSPM.XYZmm - xyz*Q) <= xY.spec(:)*Q/2));

        case 'cluster'
            %---------------------------------------------------------------
            [x i] = spm_XYZreg('NearestXYZ',xyz,xSPM.XYZmm);
            A     = spm_clusters(xSPM.XYZ);
            Q     = find(A == A(i));
    end

    % voxels defined
    %-----------------------------------------------------------------------
    spm('Pointer','Watch')

    %-Extract required data from results files
    %=======================================================================

    %-Get raw data, whiten and filter
    %-----------------------------------------------------------------------
    y        = spm_get_data(SPM.xY.VY,xSPM.XYZ(:,Q));
    y        = spm_filter(SPM.xX.K,SPM.xX.W*y);
    xY.XYZmm = xSPM.XYZmm(:,Q);


    %-Computation
    %=======================================================================

    % remove null space of contrast
    %-----------------------------------------------------------------------
    if xY.Ic

        %-Parameter estimates: beta = xX.pKX*xX.K*y
        %---------------------------------------------------------------
        beta  = spm_get_data(SPM.Vbeta,xSPM.XYZ(:,Q));

        %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
        %---------------------------------------------------------------
        y     = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);

    end

    % confounds
    %-----------------------------------------------------------------------
    xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);

    % extract session-specific rows from data and confounds
    %-----------------------------------------------------------------------
    try
        i     = SPM.Sess(xY.Sess).row;
        y     = y(i,:);
        xY.X0 = xY.X0(i,:);
    end

    % and add session-specific filter confounds
    %-----------------------------------------------------------------------
    try
        xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
    end

    %=======================================================================
    try
        xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH]; % Compatibility check
    end
    %=======================================================================

    % Remove zero columns from X0 if there are any
    % (eg. SPM.xX.iB columns may be all-zero for a particular session)
    ncol_X0=size(xY.X0,2);
    nrow_X0=size(xY.X0,1);
    new_X0=[];
    for col_X0=1:ncol_X0,
        if ~(length(find(xY.X0(:,col_X0)==0))==nrow_X0)
            new_X0=[new_X0 xY.X0(:,col_X0)];
        end
    end
    xY.X0=new_X0;


    % compute regional response in terms of first eigenvariate
    %-----------------------------------------------------------------------
    [m n]   = size(y);
    if m > n
        [v s v] = svd(spm_atranspa(y));
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
    else
        [u s u] = svd(spm_atranspa(y'));
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
    end
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    Y       = u*sqrt(s(1)/n);

    % set in structure
    %-----------------------------------------------------------------------
    xY.y    = y;
    xY.u    = Y;
    xY.v    = v;
    xY.s    = s;

    % %-Display VOI weighting and eigenvariate
    % %========================================================================
    %
    % % show position
    % %------------------------------------------------------------------------
    % spm_results_ui('Clear',Fgraph);
    % figure(Fgraph);
    % subplot(2,2,3)
    % spm_dcm_display(xY,[],[],[[1 0 0];[0 1 0]]',64)
    %
    %
    % % show dynamics
    % %------------------------------------------------------------------------
    % subplot(2,2,4)
    % try
    % 	plot(SPM.xY.RT*[1:length(xY.u)],Y)
    % 	str = 'time (seconds}';
    % catch
    % 	plot(Y)
    % 	str = 'scan';
    % end
    % title(['1st eigenvariate: ' xY.name],'FontSize',10)
    % str = {	str;' ';sprintf(...
    % 		'%d voxels in VOI at [%3.0f %3.0f %3.0f]',...
    % 		length(Q),xyz);sprintf('Variance: %0.2f%%',s(1)*100/sum(s))};
    % xlabel(str)
    % axis tight square
    %

    % save
    %-----------------------------------------------------------------------
    str     = ['VOI_' Sessname xY.name];
    if isfield(xY,'Sess')
        if length(xY.Sess) == 1
            str = sprintf('VOI_%s_%i',xY.name,xY.Sess);
        end
    end
    save(fullfile(SPM.swd,str),'Y','xY')
end

% %-Reset title
% %-----------------------------------------------------------------------
% spm('FigName',header);
% spm('Pointer','Arrow')

