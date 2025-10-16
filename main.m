%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Created by: Nili Abtahi
%
%   Laboratory for Computational Sensing and Robotics, John Hopkins University
%
%   Contact: Nili Abtahi (nabtahi1@jhu.edu)
% **************************************************************

% %
% % The main program.
% %
% % This code has five major parts:
% %    1. set geometric/physical data             ( can be changed by user     )
% %    2. prepare geometry and mesh for analysis  ( done by code automatically )
% %    3. prepare physical matrices for FEA       ( done by code automatically )
% %    4. perform modal and transient analysis    ( done by code automatically )
% %    5. post-processing data                    ( done by code automatically )
% %
% % The user can change following input data:
% %
% % 1.1. set geometric dimension of the model ( in milimeter )
% %      newLength  = 30 ;    % x-direction size : (value is 30 mm in beam_1p5x30.stl )
% %      newWidth   = 12 ;    % y-direction size : (value is 12 mm in beam_1p5x30.stl )
% %      newHeight  = 3  ;    % z-direction size : (value is 3  mm in beam_1p5x30.stl  )
% %
% % 1.2. define material property of the piezo material ( may be changed by user)
% %      E   =  2.0e9      ;  % Elastic module, N/m^2 ( default is  E   = 2.0e9     )
% %      nu  =  0.29       ;  % Poisson's ratio       ( default is  nu  = 0.29      )
% %      rho =  1.770e3    ;  % mass density, kg/ m^3 ( default is  rho = 1.770e3   )
% %
% % 1.3. Load: specified by a vector [F_x , F_y , F_z ] at tip position in free end
% %      Load = [ 0.0 , 0.0 , -2000.0 ] ;
% %
% % 1.4. Number of desired first natural frequencies to be found.
% %      numFirstNaturalFrequencies = 20 ;
% %
% %
% % 1.5. Time step and max analysis time in transient analysis.
% %      Let w1 = 2 pi /T be first natural frequency, then T = 2 * pi /w1
% %      is period.
% %      The max time for analysis is another multiple of T, i.e.  T_max = 10 * T ;
% %      The time that load acts is also a multiple of T, for example   loadTime = 0.02 * T ;
% %      The number of timeSteps is an integer, i.e. 500 ;
% % 
% %  maxTimeRatio   = 10.0    ;    % TimeInterval = [0 , maxTimeRatio * T ]; where T is fundamental period.
% %  loadTimeRatio  = 0.02    ;    % The load is acted at interval LoadTimeInterval = [0 , loadTimeRatio * T ] ;
% %  maxTimeStepRatio = 0.05  ;    % The time step begins from dt = loadTimeRatio * T / 50 ;
% %                                % and, after 80 iteration, is multiplied by 2 after every 30 iterations.
% %                                % This increment continuous as long as dt <= maxTimeStepRatio * T ;
% % % %
% % 1.6. Rayleigh damping coefficients alpha and beta,
% %      i.e. DampingMatrix = alpha * M + beta * K where M and K are mass and stiffness matrices.
% %
% %      Given two first natural frequencies w0 and w1, these parameters are
% %            alpha = zeta * ( 2 * w0 * w1 )/( w0 + w1 ) ,
% %            beta  = zeta * 2 / ( w0 + w1 )
% %      where zeta is
% %
% %      zeta  = 1/ ( 2 * Q ) ;
% %
% %      with Q being pvdf quality factor. It is given in product
% %      specification of the pvdf material. Using some articles, mentioned
% %      in document file, its value is adopted as,
% %
% %
% %      pvdfQualityFactor = 12.0 ;
% %
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars ;
clear all ;
clc ;
close all ;


tic

%% 1. Set geometric and physical parameters

% % 1.1. set geometric dimension of the model ( in milimeter )
newLength  = 30 ;    % x-direction size : (value is 30 in beam_1p5x30.stl )
newWidth   = 12 ;    % y-direction size : (value is 12 in beam_1p5x30.stl )
newHeight  = 3  ;    % z-direction size : (value is 3  in beam_1p5x30.stl  )



% % 1.2. define material property of the piezo material ( may be changed by user)
E   =  2.0e9      ;  % Elastic module, N/m^2 ( default is  E   = 2.0e9     )
nu  =  0.29       ;  % Poisson's ratio       ( default is  nu  = 0.29      )
rho =  1.770e3    ;  % mass density, kg/ m^3 ( default is  rho = 1.770e3   )

% % test case
% E   =  1.56e9      ;  % Elastic module, N/m^2 
% nu  =  0.20        ;  % Poisson's ratio     
% rho =  1.96211e3   ;  % mass density, kg/ m^3 


% % 1.3. Load: specified by a vector [F_x , F_y , F_z ] at tip position in free end ( may be changed by user)
Load = [ 0.0 , 0.0 , -100.0 ] ;


% % 1.4. Number of desired first natural frequencies to be found.
numFirstNaturalFrequencies = 10 ;


% % 1.5. Time step and max analysis time in transient analysis.
% %      Let w1 = 2 pi /T be first natural frequency, then T = 2 * pi /w1
% %      is period.
% %      The max time for analysis is another multiple of T, i.e.  T_max = 10 * T ;
% %      The time that load acts is also a multiple of T, for example   loadTime = 0.02 * T ;
% %      The number of timeSteps is an integer, i.e. 500 ;

 maxTimeRatio   = 10.0   ;    % TimeInterval = [0 , maxTimeRatio * T ]; where T is fundamental period.
 loadTimeRatio  = 0.02   ;    % The load is acted at interval LoadTimeInterval = [0 , loadTimeRatio * T ] ;
 maxTimeStepRatio = 0.03 ;    % The time step begins from dt = loadTimeRatio * T / 50 ;
                              % and, after 80 iteration, is multiplied by 2 after every 30 iterations
                              % as long as dt <= maxTimeStepRatio * T ;
 
% % 1.6. pvdf quality factor

 pvdfQualityFactor = 12.0 ;


 %% 2. Prepare geometry and physics of the model

    geometry = mGeometry("beam_5x30.stl") ;                           %% 2.1. read geometry
    geometry.rescaleGeometry([ newLength , newWidth , newHeight ] ) ; %% 2.2. rescale geometry
    geometry.specifyTipVertex() ;                                     %% 2.3. specify tip vertex
    geometry.createMesh()  ;                                          %% 2.4. create mesh
    geometry.fixFaceAt('xmax') ;      %% 2.5. fix the plane (face) specified by x = xmax is fixed

  
    
  %% 3. Prepare physics

   physics = mPhysics(geometry) ;    %% 3.1. create physics object

   physics.setDofMatrix(    geometry ) ;   %% 3.2. set dof matrix


   physics.assmbleMatrices( geometry , E , nu , rho , Load ) ;   %% 3.3. do FE analysis to form global mass and stiffness matrices, and force vector


   physics.freeUnusedMemory() ;

   physics.doModalAnalysis( numFirstNaturalFrequencies ) ;

%%
   zeta  = 1.0 / ( 2.0 * pvdfQualityFactor ) ;


   alpha =   zeta * 2.0 * physics.naturalFrequency(1) * physics.naturalFrequency(2) /...
                        ( physics.naturalFrequency(1) + physics.naturalFrequency(2) ) ;

   beta0   = ( zeta * 2.0 ) / ( physics.naturalFrequency(1) + physics.naturalFrequency(2) ) ;

   physics.doTransientAnalysis( alpha , beta0 ,  maxTimeRatio , loadTimeRatio , maxTimeStepRatio ) ;

   %%
   physics.fitToDecayingCurve( geometry , Load  ) ;

   physics.destrucPhysics() ;
   geometry.destrcutGeometry() ;
   toc



