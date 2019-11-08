function [Gfunct_r, Gfunct_theta,Gfunct_phi] = calcGreensFunction_weighted_correct_areas( sc_pos_r,sc_pos_theta,sc_pos_phi, rsurface,thetas_surf, phis_surf,arealist)
%Calculate the Green's function to upward continue 
%the surface Br field at surface points: (rsurface, thetas_surf,phis_surf) 

%to satellite points: sc_pos_r, sc_pos_theta,sc_pos_phi.

% where:
% rsurface is constant (3393.5 km for the Mars models in this paper)
% thetas & phis are in radians; phi is E longitude and theta is co-latitude
% arealist should sum to 4pi.

%Br_satellite     = Gfunct_r     *Br_surf
%Btheta_satellite = Gfunct_theta *Br_surf
%Bphi_satellite   = Gfunct_phi   *Br_surf


%declare/allocate constants
fourpi = 4*pi;
c = rsurface; %assumes that we are using a constant radius of the surface and not an oblate planet.
Gfunct_r = zeros(length(sc_pos_theta), length(thetas_surf));
Gfunct_theta = zeros(length(sc_pos_theta), length(thetas_surf));
Gfunct_phi = zeros(length(sc_pos_theta), length(thetas_surf));


for zz = 1:length(sc_pos_theta) %loop over satellite positions
    a = sc_pos_r(zz);
    b = c/a;
    b2 =  b^2;
    b3 = b^3;
    const1  = b2*(1-b2)/fourpi;
    constt2 = b3/fourpi;
    
    %pre-compute things to save time
    the_sat = sc_pos_theta(zz);
    ph_sat = sc_pos_phi(zz);
    costh_sat = cos(the_sat);
    sinth_sat = sin(the_sat);
    cosph_sat = cos(ph_sat);
    sinph_sat = sin(ph_sat);
    
    %loop over surface positions
    for ii = 1:length(thetas_surf)
        the_surf = thetas_surf(ii);
        costh_surf = cos(the_surf);
        sinth_surf = sin(the_surf);
        costh_sat_times_costh_surf = costh_sat*costh_surf;
        sinth_sat_times_sinth_surf = sinth_sat*sinth_surf;
        costh_sat_times_sinth_surf = costh_sat*sinth_surf;
        sinth_sat_times_costh_surf = sinth_sat*costh_surf;
        
        ph_surf = phis_surf(ii);
        cosph_surf = cos(ph_surf);
        sinph_surf = sin(ph_surf);
        
        mu = costh_sat_times_costh_surf + sinth_sat_times_sinth_surf*cos(ph_sat-ph_surf); 
        dmu_dphi   = sinth_sat_times_sinth_surf*( -cosph_surf *sinph_sat + sinph_surf*cosph_sat   );
        dmu_dtheta = sinth_sat_times_costh_surf  -costh_sat_times_sinth_surf*    (cosph_surf*cosph_sat + sinph_surf*sinph_sat)    ;
        
        f = (1-2*b*mu+b2)^0.5;
        Tt = 1+f-mu*b;
        
        
        %save into Gfunct
        Gfunct_r(zz,ii) =  const1/(f^3)*arealist(ii);  %This is the  Gubbins way from the paper, for Gr
        Gfunct_phi(zz,ii) = -(constt2*(1+2*f-b2)/(f^3*Tt))*dmu_dphi/sinth_sat*arealist(ii);   
        Gfunct_theta(zz,ii) = (constt2*(1+2*f-b2)/(f^3*Tt))*dmu_dtheta*arealist(ii);
    end
end



end

