clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Moore et al.,2017, JGR Planets
%%%%%%%% This code loads a model of the field from my paper, and plots it
%%%%%%%% in various ways. It also explains all the variable names and how
%%%%%%%% to use them.
%%%%%%%% please cite Moore et al., 2017, JGR planets when acknowledging use
%%%%%%%% of these models or code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load coordinates:
load('jgre20703-sup-0002-supinfo.mat')
%%%%%%% variables: %%%%%%%%%%%%%%%%%%
% this file contains the coordinates of sphere which has been tesselated
% into 10,000 cells using a voronoi method
%     voronoi boundary: this list contains the boundary points defining each of
%            the 10,000 cells
%     V: x,y,z coordinates for the center of the 10,000 cells (defined on a unit sphere of radius = 1)
%     arealist: the area of each voronoi cell. relevant for calculation of the
%             field using a Green's function method
rMars = 3393.5; %km. 
xs = rMars*V(:,1);
ys = rMars*V(:,2);
zs = rMars*V(:,3);

rs_surf     = rMars* sqrt(xs.^2+ys.^2+ zs.^2);
thetas_surf = atan2((xs.^2+ys.^2).^0.5,zs);
phis_surf   = atan2(ys,xs);
for ii = 1:length(phis_surf)
    if(phis_surf(ii)<0)
        phis_surf(ii) = phis_surf(ii)+2*pi;
    end
end
    
clearvars V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load magnetic field model
%load('jgre20703-sup-0004-supinfo.mat') %%% L1
%load('jgre20703-sup-0005-supinfo.mat') %%% L1
%load('jgre20703-sup-0006-supinfo.mat') %%% L1
load('jgre20703-sup-0007-supinfo.mat') %%% elastic net
%%%%%%% variables: %%%%%%%%%%%%%%%%%%
% these files each contain a different magnetic field model calculated using
% different values of the regularization parameters
%     alpha: controls the amount of L1 vs L2 norm used in the inversion
%         (ranges from 0 to 1. alpha = 1 means a pure L1 norm)
%     lamba1: controls the relative weighting of the norm vs the misfit in
%         the inversion
% the resulting model, B_surf_glmnet, is a list of the magnetic field
% values for each of the 10,000 grid cells of a sphere (as defined in the
% coordinates above). It is assumed that the field is constant over the
% entirety of a cell. 
% the model has the resulting properties of
%     percent_zero: percent of grid cells with a surface Br field of 0.
%     rms_misfit_nT_glmnet: rms misfit (in nT) to the original MGS dataset
%     Br_sat_glmnet: the radial magnetic field predicted by my model at the coordinates of the original
%           MGS dataset (in nT)
Btheta_sat_glmnet = B_sat_glmnet(80230+1  :2*80230);
Bphi_sat_glmnet   = B_sat_glmnet(2*80230+1:3*80230);
%     Btheta_sat_glmnet: the theta-component of the magnetic field predicted by my model at the coordinates of the original
%         MGS  (in nT)
%     Bphi_sat_glmnet: the phi-component of the magnetic field predicted by my model at the coordinates of the original
%         MGS dataset (in nT)
%         ***All coordinates and fields listed here are defined using a right-handed
%         coordinate system (i.e. East longitudes instead of West
%         longitudes). ***
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the magnetic field model at Mars' surface
disp('1. Plotting surface magnetic field from model...(Figure 1)')
figure(1)
subplot(2,1,1)
scatter3(xs,ys,zs,10,Br_surf_glmnet,'filled')
colorbar
colormap('jet')
title('Surface Br, in nT (Moore et al., 2017)')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(2,1,2)
scatter(phis_surf*180/pi, thetas_surf*180/pi, 5, Br_surf_glmnet, 'filled')
set(gca,'YDir','reverse')
axis([0 360 0 180])
%caxis([-150 325])
colorbar
colormap('jet')
title('Moore et al. model, surface Br(nT)')
xlabel 'E Longitude (deg)'; ylabel 'Co-latitude (deg)'
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load original MGS satellite data 
load('jgre20703-sup-0003-supinfo.mat')
%%%%%%% variables: %%%%%%%%%%%%%%%%%%
% This file contains the MGS dataset I used in my inversion. For details on
% data selection procedures, please see Moore et al., 2017
%    sc_pos_r: Radial coordinate of the MGS spacecraft (distance from
%          center of Mars), in km
%    sc_pos_theta: theta-component of the MGS spacecraft's position, in
%          radians
%    sc_pos_phi: phi-component of the MGS spacecraft's position, in radians
%         ***All coordinates and fields listed here are defined using a right-handed
%         coordinate system (i.e. East longitudes instead of West
%         longitudes). ***
%    gamma_br: radial magnetic field (nT) at each location
%    gamma_btheta: theta-component of the magnetic field (nT) at each location
%    gamma_bphi: phi-component of the magnetic field (nT) at each location
%    year: corresponding year of each datapoint
%    decimal_day: day of year
%    sao_i: a parameter describing the state of the solar panel. used to find "night-side" data. see PDS for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert satellite coordinates to x,y,z:
sc_pos_x = sc_pos_r.*sin(sc_pos_theta).*cos(sc_pos_phi);
sc_pos_y = sc_pos_r.*sin(sc_pos_theta).*sin(sc_pos_phi);
sc_pos_z = sc_pos_r.*cos(sc_pos_theta);

disp('2. Plotting surface magnetic field from model...(Figure 2)')
figure(2)
%%%% plot satellite field at MGS locations
subplot(2,1,1)
scatter(sc_pos_phi*180/pi, sc_pos_theta*180/pi, 20, gamma_br, 'filled')
set(gca,'YDir','reverse')
axis([0 360 0 180])
caxis([-150 325])
colorbar
colormap('jet')
title('MGS data, Br at satellite altitude (nT)')
xlabel 'E Longitude (deg)'; ylabel 'Co-latitude (deg)'
drawnow

%%%%% plot Moore et al. field model at MGS locations
subplot(2,1,2)
scatter(sc_pos_phi*180/pi, sc_pos_theta*180/pi, 20, Br_sat_glmnet, 'filled')
set(gca,'YDir','reverse')
axis([0 360 0 180])
caxis([-150 325])
colorbar
colormap('jet')
title('Moore et al. model, Br at satellite altitude (nT)')
xlabel 'E Longitude (deg)'; ylabel 'Co-latitude (deg)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you'd like to calculate the predicted field at a different location using any of the Moore at
% al. models, you can use the Green's function function:
disp('3. Using a Greens function to calculate the field at a new altitude')

%%%% example coordinates you'd like to plot the model at:
desired_radius = rMars + 200; %e.g. what does the field look like at 200 km??
npoints        = 5000; % let's use 1000 random points at the desired radius (200km altitude)
rcalc          = desired_radius*ones(npoints,1);
thetas_calc    = pi  *rand(npoints,1); %colatitudes, in radians
phis_calc      = 2*pi*rand(npoints,1);

%%%%% just put your new coordinates into the Green's function (included).
%%%%%     Thetas_surf and phis_surf are the coordinates of your surface grid,
%%%%%     arealist is the surface area of each grid cell (in terms of a unit sphere)
%%%%%          (from the very first file we loaded).
%%% Be careful how many points you use in npoints, as these matrices can
%%% get quite large. 
disp('......calculating Greens function...')
size_of_each_Gfunct_component_matrix_in_GB = npoints*1000*8/1e9
tic
[Gfunct_r, Gfunct_theta,Gfunct_phi] = calcGreensFunction_weighted_correct_areas( rcalc,...
                                         thetas_calc,phis_calc , rMars,thetas_surf, phis_surf,arealist);
toc
%%%%% Now use the Green's function to calculate the field at altitude:
Br_sat_desired_radius     = Gfunct_r    *Br_surf_glmnet;
Btheta_sat_desired_radius = Gfunct_theta*Br_surf_glmnet;
Bphi_sat_desired_radius   = Gfunct_phi  *Br_surf_glmnet;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot your results:
disp('......plotting (Figure 3)')
figure(3)
subplot(3,1,1)
scatter(phis_calc*180/pi, thetas_calc*180/pi, 10, Br_sat_desired_radius , 'filled')
set(gca,'YDir','reverse')
axis([0 360 0 180])
colorbar
colormap('jet')
title('Moore et al. model, Br at new points (nT)')
xlabel 'E Longitude (deg)'; ylabel 'Co-latitude (deg)'

subplot(3,1,2)
scatter(phis_calc*180/pi, thetas_calc*180/pi, 10, Btheta_sat_desired_radius , 'filled')
set(gca,'YDir','reverse')
axis([0 360 0 180])
colorbar
colormap('jet')
title('Moore et al. model, Btheta at new points (nT)')
xlabel 'E Longitude (deg)'; ylabel 'Co-latitude (deg)'

subplot(3,1,3)
scatter(phis_calc*180/pi, thetas_calc*180/pi, 10, Bphi_sat_desired_radius , 'filled')
set(gca,'YDir','reverse')
axis([0 360 0 180])
colorbar
colormap('jet')
title('Moore et al. model, Bphi at new points (nT)')
xlabel 'E Longitude (deg)'; ylabel 'Co-latitude (deg)'


