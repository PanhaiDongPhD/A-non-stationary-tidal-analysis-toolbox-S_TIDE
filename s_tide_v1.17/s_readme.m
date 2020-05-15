% s_readme.m
%
%                    Tidal Analysis Toolbox S_TIDE
%                              by        
%               Haidong Pan (Ocean University of China) 
%
%                 Version 1.17    May 10, 2020
%
%Please citing:
%1. Pan, H., X. Lv, Y. Wang, P. Matte, H. Chen, and G. Jin (2018), Exploration of Tidal-Fluvial Interaction in the Columbia River Estuary Using S_TIDE, J. Geophys. Res. Ocean., 123(9), 6598–6619, doi:10.1029/2018JC014146.
%2.	Jin, G., H. Pan, Q. Zhang, X. Lv, W. Zhao, and Y. Gao (2018), Determination of Harmonic Parameters with Temporal Variations: An Enhanced Harmonic Analysis Algorithm and Application to Internal Tidal Currents in the South China Sea, J. Atmos. Ocean. Technol., 35(7), 1375-1398
%3. Pawlowicz, R., B. Beardsley, and S. Lentz (2002), Classical tidal harmonic analysis including error estimates in MATLAB using T_TIDE, Comput. Geosci., 28(8), 929–937
%
% The toolbox presently contains the following mfiles:
%
% ---FOR ANALYSIS
% s_tide.m       - computes the tidal analysis of the real tidal 
%                  time series using Enhanced Harmonic Analysis (EHA).
% s_nodal_correction.m   -correct the amplitues and phases.
%
% t_vuf.m        - computes nodal corrections (from T_TIDE package). 
% 
% t_astron.m     - computes astronomical arguments (from T_TIDE package). 
%
% t_getconsts.m     - Gets constituent data structures (from T_TIDE package).
%
% s_calculate_coefficient  -computes the spline interpolation weights
%
% l_calculate_coefficient  -computes the linear interpolation weights
%
% sinc_calculate_coefficient  -computes the sinc interpolation weights
%   
% s_tide_m1.m                 - a modified version of s_tide.m
%
% s_tide_m2.m                 - a modified version of s_tide.m
%
% s_tide_m3.m                 - a modified version of s_tide.m (for satellite data)
%
% ---FOR DOCUMENTATION
%
% s_readme.m     - this file.
%
% ---FOR DEMONSTRATION
%
% s_demo.m       - a short example using the Kushiro elevation data.
% S_TIDE toolbox tutorial/S_TIDE工具包中文教程
% 
%
%
% Various data files are also included:
%
% t_constituents.mat - constituent data structures.
%
% kushiro.mat  -  hourly elevation data set used as an 
%                  example (obtained from PSMSL)
% imf.mat - the modes obtained by Empirical Mode Decomposition (EMD) of the
%           kushiro elevation data.
% K1_amp_hilo.mat - monthly K1 tidal amplitudes at Hilo
% N2_amp_astoria.mat  - monthly N2 tidal amplitudes at Astoria
% wt_SIO.mat - daily water temperature at Sion
% satellite.mat - sea levels at 14.1633N,113.6244E obtained from satellite(T/P,J1,J2)
%
% Questions or comments to:
% Haidong Pan(Ocean University of China) 2017/05/02   潘海东（中国海洋大学）
% email:1027128116@qq.com    panhaidong_phd@qq.com
% S_TIDE users QQ group number (QQ群号）485997351

help s_readme