clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Moore et al.,2017, JGR Planets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load model and plot at surface:

%%%%%%%% load coordinates:
load('jgre20703-sup-0002-supinfo.mat')
rMars = 3393.5;
xs = rMars*V(:,1);
ys = rMars*V(:,2);
zs = rMars*V(:,3);


%%%%%%%% load model and plot:
%load('jgre20703-sup-0004-supinfo.mat') %%% L1
%load('jgre20703-sup-0005-supinfo.mat') %%% L1
%load('jgre20703-sup-0006-supinfo.mat') %%% L1
load('jgre20703-sup-0007-supinfo.mat') %%% elastic net

figure(1)
scatter3(xs,ys,zs,20,Br_surf_glmnet,'filled')
colorbar
colormap('jet')
title('Br at the planetary surface, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load satellite data and plot model against data
load('jgre20703-sup-0003-supinfo.mat')
figure(2)

%convert satellite coordinates to x,y,z:
sc_pos_x = sc_pos_r.*sin(sc_pos_theta).*cos(sc_pos_phi);
sc_pos_y = sc_pos_r.*sin(sc_pos_theta).*sin(sc_pos_phi);
sc_pos_z = sc_pos_r.*cos(sc_pos_theta);


%%%%% plot satellite dataset:
subplot(3,3,1)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,gamma_br,'filled')
colorbar
colormap('jet')
title('MGS data, Br at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(3,3,2)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,gamma_btheta,'filled')
colorbar
colormap('jet')
title('MGS data, Btheta at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(3,3,3)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,gamma_bphi,'filled')
colorbar
colormap('jet')
title('MGS data, Bphi at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';


%%%%% plot model:
len = length(Br_sat_glmnet);
Btheta_sat_glmnet = B_sat_glmnet(len+1  :2*len);
Bphi_sat_glmnet   = B_sat_glmnet(2*len+1:3*len);

subplot(3,3,4)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,Br_sat_glmnet,'filled')
colorbar
colormap('jet')
title('Model, Br at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(3,3,5)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,Btheta_sat_glmnet,'filled')
colorbar
colormap('jet')
title('Model, Btheta at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(3,3,6)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,Bphi_sat_glmnet,'filled')
colorbar
colormap('jet')
title('Model, Btheta at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

%%%%% plot residuals:
subplot(3,3,7)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,gamma_br-Br_sat_glmnet,'filled')
colorbar
colormap('jet')
title('Residual, Br at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(3,3,8)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,gamma_btheta-Btheta_sat_glmnet,'filled')
colorbar
colormap('jet')
title('Residual, Btheta at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';

subplot(3,3,9)
scatter3(sc_pos_x,sc_pos_y,sc_pos_z,20,gamma_bphi-Bphi_sat_glmnet,'filled')
colorbar
colormap('jet')
title('Residual, Btheta at satellite altitude, in nT')
xlabel 'x in km'; ylabel 'y in km'; zlabel 'z in km';



