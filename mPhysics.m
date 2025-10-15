%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  The mPhysics does following jobs:
% %     1. forming degree of freedom matrix,
% %     2. calculating global mass, stiffness matrices and force vector,
% %     3. doing modal analysis for calculating natural frequencies,
% %     4. doing a transient analysis in time domain.
% %
% %  This class uses only basic matlab functions, i.e. no FEM-based
% %  function in PDE toolbox is used in this class.
% %  This code needs Matlab R2015+
% %
% % The displacement/velocity/acceleration array stores data at given nodes.
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef  mPhysics < handle

    %%==============================================================================%%
     properties
        displacement        = []  ;  % displacement array at given nodes (tip and others) for all time steps
        naturalFrequency    = []  ;  % natural frequency in cycle / sec

     end  % of properties


    properties (Access = private)
        numVars             = 0   ;  % number of variables
        mass                = []  ;  % global mass matrix
        stiffness           = []  ;  % global stiffness matrix
        force               = []  ;  % global force vector
        degreeOfFreedom     = []  ;  % degree of freedom matrix
        elementalMass       = []  ;  % elemental mass matrix
        elementalStiffness  = []  ;  % elemental stiffness matrix
        rowVarSet           = []  ;  % row-to-row local to global index convertion
        shapeFunctionVector = []  ;  % vector containing shape function
        dLocalShapeFuns     = []  ;  % derivative of shape functions w.r.t. local variables
        JacobianMatrix      = []  ;  % Jacobain matrix of global to local transformation
        inverseJacobianT    = []  ;  % inverse of transpose of Jacobain matrix
        dGlobalShapeFuncs   = []  ;  % derivative of shape functions w.r.t. global variables.
        rigidity            = []  ;  % Rigidity matrix of material
        BmatrixMu           = []  ;  % matrix used in elemental stiffness matrix
        BmatrixNu           = []  ;  % matrix used in elemental stiffness matrix
        tipNodeId           = 0   ;  % index of nearesr node to tip-vertex
        permutation         = []  ;  % permutation matrix
        equivStiffness      = []  ;  % equivalent stiffness matrix
        equivForce          = []  ;  % equivalent force
        equivSolution       = []  ;  % solution of equivalent system.
        tempVec             = []  ;  % a temporary vector
        time                = []  ;
        
        % % 15 Gauss nodes for Gauss quadrature
        gaussPoint = [ 0.2500000000000000, 0.2500000000000000, 0.2500000000000000 ;
                       0.0000000000000000, 0.3333333333333333, 0.3333333333333333 ;
                       0.3333333333333333, 0.3333333333333333, 0.3333333333333333 ;
                       0.3333333333333333, 0.3333333333333333, 0.0000000000000000 ;
                       0.3333333333333333, 0.0000000000000000, 0.3333333333333333 ;
                       0.7272727272727273, 0.0909090909090909, 0.0909090909090909 ;
                       0.0909090909090909, 0.0909090909090909, 0.0909090909090909 ;
                       0.0909090909090909, 0.0909090909090909, 0.7272727272727273 ;
                       0.0909090909090909, 0.7272727272727273, 0.0909090909090909 ;
                       0.4334498464263357, 0.0665501535736643, 0.0665501535736643 ;
                       0.0665501535736643, 0.4334498464263357, 0.0665501535736643 ;
                       0.0665501535736643, 0.0665501535736643, 0.4334498464263357 ;
                       0.0665501535736643, 0.4334498464263357, 0.4334498464263357 ;
                       0.4334498464263357, 0.0665501535736643, 0.4334498464263357 ;
                       0.4334498464263357, 0.4334498464263357, 0.0665501535736643 ] ;


     gaussWeight  = [ 0.1817020685825351, 0.0361607142857143, 0.0361607142857143, 0.0361607142857143, 0.0361607142857143,...
     0.0698714945161738, 0.0698714945161738, 0.0698714945161738, 0.0698714945161738, 0.0656948493683187,...
     0.0656948493683187, 0.0656948493683187, 0.0656948493683187, 0.0656948493683187, 0.0656948493683187]/6.0 ;

%        gaussPoint  = [ 0.1381966011250105 , 	0.1381966011250105 , 0.1381966011250105 ;
%                        0.5854101966249685 ,  0.1381966011250105 , 0.1381966011250105 ;
%                        0.1381966011250105 ,  0.5854101966249685 , 0.1381966011250105 ;
%                        0.1381966011250105 ,  0.1381966011250105 , 0.5854101966249685 ] ; % 4-point Gauss integration
%        gaussWeight = [ 1.0 /24.0 , 1.0 /24.0 , 1.0 /24.0 , 1.0 /24.0     ] ; % weight of each term in integration
        fitParameters       = [] ;
    end


    % % public methods
    methods

          %%=======================================================================================%%
          function obj = mPhysics( geometryObj )  %% constructor
              if ( nargin == 1 )
                  if( isa(geometryObj ,'mGeometry') )
                      obj.degreeOfFreedom = ones( size( geometryObj.nodes ) , 'int32') ;
                  else
                      error("The constructor of physics class needs a geometry object as input")
                  end
              else
                  error("The constructor of physics class needs a geometry object as input")
              end
          end % of constructor

          %%=======================================================================================%%
          function setDofMatrix( obj , geomObj )  % % set degree of freedom matrix
                id  = geomObj.fix_face_axis ;
                val = geomObj.fix_face_value ;
                minLength = geomObj.minElementSize ;
                obj.degreeOfFreedom( abs( geomObj.nodes(:,id) - val ) < 1.e-3 * minLength , : ) = 0 ;

                obj.numVars = 0 ;
                for i = 1 : geomObj.numNodes                         % % loop over num nodes
                    for j = 1:3                                      % % loop over numDim vars = [u_x , u_y, u_z]
                        if( obj.degreeOfFreedom(i,j) == 1 )          % % check if node has dof (i.e. not fixed)
                            obj.numVars = obj.numVars + 1 ;          % % a new variable is found
                            obj.degreeOfFreedom(i,j) = obj.numVars ; % % that variable is labeled by its global index, i.e. numVars
                        end
                    end % of for j-loop
                end     % of for i-loop

          end % of setDofMatrix function

          %%=======================================================================================%%
          function  assmbleMatrices( obj , geomObj , E , nu , rho , Load )     % % calculate global mass and stiffness matrix


               % % 1. allocate required matrices.
               % % 1.1. assign global matrices

               obj.mass      = zeros( [ obj.numVars , obj.numVars ] , 'double'  ) ;    % % square matrix with numVars cols and rows
               obj.stiffness = zeros( [ obj.numVars , obj.numVars ] , 'double'  ) ;    % % square matrix with numVars cols and rows
               obj.force     = zeros( [obj.numVars, 1 ] , 'double'  ) ;                % %  matrix with numVars rows and 3 cols

               % % 1.2. assign elemental matrices. The elements are quadratic tetrahedral having 10 nodes, each of which 3 dof, so 30 vars.
               obj.elementalMass        = zeros([30 , 30] , 'double' ) ;
               obj.elementalStiffness   = zeros([30 , 30] , 'double' ) ;
               obj.rowVarSet            = zeros([30 , 1 ] , 'int32'  ) ;
               obj.shapeFunctionVector  = zeros([1  , 10] , 'double' ) ;
               obj.dLocalShapeFuns      = zeros([3  , 10] , 'double' ) ;
               obj.JacobianMatrix       = zeros([3  , 3 ] , 'double' ) ;
               obj.inverseJacobianT     = zeros([3  , 3 ] , 'double' ) ;
               obj.rigidity             = zeros([6  , 6 ] , 'double' ) ;
               obj.BmatrixMu            = zeros([6  , 3 ] , 'double' ) ;
               obj.BmatrixNu            = zeros([6  , 3 ] , 'double' ) ;

               % % 1.3. set Rigidity matrix as it is constant for all elements
               t0    = E / ( ( 1.0 - 2.0 * nu ) * ( 1.0 + nu ) ) ;
               obj.rigidity(1,1)   = t0 * ( 1.0 - nu )  ;
               obj.rigidity(2,2)   = obj.rigidity(1,1) ;
               obj.rigidity(3,3)   = obj.rigidity(1,1) ;
               obj.rigidity(1,2)   = t0 * nu ;
               obj.rigidity(1,3)   = obj.rigidity(1,2)  ;
               obj.rigidity(2,1)   = obj.rigidity(1,2)  ;
               obj.rigidity(2,3)   = obj.rigidity(1,2)  ;
               obj.rigidity(3,1:2) = obj.rigidity(1,2)  ;
               obj.rigidity(4,4)   = t0 * ( 0.5 - nu )  ;
               obj.rigidity(5,5)   = obj.rigidity(4,4)  ;
               obj.rigidity(6,6)   = obj.rigidity(4,4)  ;

               % % 1.4. find nearest node to tip vertex
               minDist   = inf ;
               obj.tipNodeId = 0 ;

               for i = 1 : geomObj.numNodes
                   dr = geomObj.tip_position - geomObj.nodes(i,:) ;
                   dist = norm( dr) ;
                   if( dist < minDist )
                       minDist = dist ;
                       obj.tipNodeId = i ;
                   end
               end % if for-i loop

               % % set force vector
               for i = 1 : 3
                   obj.force( obj.degreeOfFreedom( obj.tipNodeId , i )  ) = Load(i) ;
               end


               % % Now change geometry coordinates to SI unit, i.e. they are primarily in milimeter, change to meter
               geomObj.nodes = geomObj.nodes / 1000.0 ;  % % change from mm to meter (SI unit of length)

               % % 1.5. Do FEM for each element

               for i = 1 : geomObj.numElements
                   disp(strcat( ' ** Assemblage: element id =  ' , num2str(i) ) )
                   obj.doFEMforElement( geomObj , i , rho ) ;
               end % of for i-loop


          end % of assmbleMatrices


          %%=======================================================================================%%
          function freeUnusedMemory(obj)
               obj.elementalMass       = []  ;  % elemental mass matrix
               obj.elementalStiffness  = []  ;  % elemental stiffness matrix
               obj.rowVarSet           = []  ;  % row-to-row local to global index convertion
               obj.shapeFunctionVector = []  ;  % vector containing shape function
               obj.dLocalShapeFuns     = []  ;  % derivative of shape functions w.r.t. local variables
               obj.JacobianMatrix      = []  ;  % Jacobain matrix of global to local transformation
               obj.inverseJacobianT    = []  ;  % inverse of transpose of Jacobian matrix
               obj.dGlobalShapeFuncs   = []  ;  % derivative of shape functions w.r.t. global variables.
               obj.rigidity            = []  ;
               obj.BmatrixMu           = []  ;
               obj.BmatrixNu           = []  ;
               obj.gaussPoint          = []  ;
               obj.gaussWeight         = []  ;
          end % of destructor function

          %%=======================================================================================%%
          function destrucPhysics(obj)
               obj.displacement      = []  ;  % displacement array at given nodes (tip and others) for all time steps
               obj.naturalFrequency  = []  ;  % acceleration array at given nodes (tip and others) for all time steps
               obj = [] ;
          end % of destructor function

          %%=======================================================================================%%
          function doModalAnalysis( obj , numFreqs  )
              disp('********** Beginning modal analysis' )

              % % 1. get eigenvalues of system: ( K * v = w^2 * M * v )  with w being angular frequency.
              [ ~ , d , flag] = eigs( obj.stiffness , obj.mass , numFreqs , 'smallestabs'  ) ;
              if( flag ~= 0 )
                  error(' Problem in finding natural frequencies') ;
              end


              % % 2. get natural frequencies: eigenvalues are w^2
              obj.naturalFrequency = zeros( [numFreqs,1] ) ;
              for i = 1:numFreqs
                  obj.naturalFrequency(i) = sqrt( abs( d(i,i) ) ) ;
              end
              % % 4. Write natural frequencies in text file

              fileID = fopen('./output/NaturalFrequencies.txt','w') ;
              fprintf(fileID,'Mode Id,   Frequency (cycle/second) \n'  ) ;
              fprintf(fileID,'\n'  );
              for i = 1 : numFreqs
                  fprintf(fileID,'%3d \t  %15.12f \n' , i , obj.naturalFrequency(i) / ( 2.0 * pi )  ) ;
              end
              fclose(fileID);

              disp('********** Ending modal analysis' )
          end % of function doModalAnalysis


          %%=====================================================================================%%
          function  doTransientAnalysis( obj , alpha , beta0 ,  maxTimeRatio , loadTimeRatio , maxTimeStepRatio )

              % % 0. Time step and time domain
              w1 = obj.naturalFrequency(1) ;
              T  = 2.0 * pi / w1 ;   % fundamental period
              
              % % 0. set time vector
              % % 0.1. Initially dt = loadTime / 50.0 for first 100 iteration 
              % %      If final time is reached then smaller number of iterations is used.
              
              N = maxTimeRatio / maxTimeStepRatio ;
              n_vec = (0:N)' ;
              delta1 = maxTimeStepRatio * T ;
              delta0 = loadTimeRatio * T / 40.0  ;
              b = ( delta1 - delta0 ) / ( N - 1 ) ;
              a = delta0 - b ;
              obj.time = n_vec .* ( a + b * n_vec ) ;
              
              loadTime    = loadTimeRatio * T ;
              
              
%               timeStep = loadTime / 40.0 ;
%               obj.time = (0.0:timeStep: min( 80.0 * timeStep ,  maxTime) )' ;
%               dt = ( maxTimeStep - timeStep ) / 370.0 ;
%               obj.time = [obj.time ; ( obj.time(end) + dt : dt : obj.time(end) + 370.0 * dt )' ] ;
              
%               % % 0.2. Then dt is multiplied by two for each 20 iterations until dt <= maxTimeStep
%               if( 80.0 * timeStep < maxTime )
%                   while( obj.time(end) < maxTime )
%                       timeStep = min( 2.0 * timeStep , maxTimeStep ) ;
%                       lowerBound = obj.time(end) ;
%                       [ upperBound , minId] = min( [ lowerBound + 30.0 * timeStep , maxTime ] ) ;
%                       obj.time = [obj.time ; ( lowerBound + timeStep : timeStep : upperBound )' ] ;
%                       if( minId == 2 )
%                           break ;
%                       end 
%                   end % of while
%               end % of if
              
              % % 1. define coefficients for constant acceleration Newmark method
              betaNewmark  = 0.25 ;
              gammaNewmark = 0.50 ;
              

              % % 2. allocate memory for storing data at tip node for different time steps

              obj.displacement   = zeros( [ length(obj.time) , 3] , 'double' ) ; %  displacement vector

              % % 3. Initial displacement and velocity are zero, define U, velocity and acceleration at two steps
              obj.equivForce     = zeros( [obj.numVars , 1]  , 'double' ) ;
              obj.equivStiffness = zeros( [obj.numVars , obj.numVars] , 'double' ) ;
              obj.equivSolution  = zeros( size( obj.equivForce ) , 'double' ) ;
              obj.tempVec        = zeros( [1 , obj.numVars]  , 'double' ) ;

              U              = zeros( size( obj.equivForce ) , 'double' ) ;
              Udot           = zeros( size( obj.equivForce ) , 'double' ) ;
              Uddot          = zeros( size( obj.equivForce ) , 'double' ) ;
              UdotNew        = zeros( size( obj.equivForce ) , 'double' ) ;
              UddotNew       = zeros( size( obj.equivForce ) , 'double' ) ;
              tempU          = zeros( size( obj.equivForce ) , 'double' ) ;



              % % 4. Do time domain analysis
              for i = 1 :  length(obj.time) -1
                  %disp( strcat( ' ** Transient analysis at step =  ' , num2str(i) ) )
                  fprintf('******* Transient analysis at step = %d (of %d steps) \n\n',i , length( obj.time) ) ;
                  % % 4.0 set timeStep and calculate newmark coeffs
                  
                  timeStep = obj.time(i+1) - obj.time(i) ;
                  alpha0 = 1.0 / ( betaNewmark  * timeStep^2 ) ;
                  alpha1 = gammaNewmark   / ( betaNewmark  * timeStep   ) ;
                  alpha2 = 1.0 / ( betaNewmark  * timeStep   ) ;
                  alpha3 = 1.0 / (2.0 * betaNewmark  ) - 1.0 ;
                  alpha4 = gammaNewmark   / betaNewmark  - 1.0 ;
                  alpha5 = ( gammaNewmark / ( 2.0 * betaNewmark ) - 1.0 ) * timeStep ;

                  fprintf('Time step = %12.8f ,  current time = %12.8f ,  final time = %12.8f  \n',...
                                timeStep , obj.time(i+1) , obj.time(end) ) ;
                            
                 % % 4.1 Effective stiffness matrix is not constant in all steps
                 obj.equivStiffness = ( 1.0 + alpha1 * beta0 ) * obj.stiffness + ( alpha0 + alpha1 * alpha ) * obj.mass  ;

              
                  % % 5.1. set effective force vector to zero
                  obj.equivForce(:) = 0.0 ;

                  % % 4.2. Form effective force
                  % % 4.2.1. If counter is smaller than numStepLoad then load is added
                  if( obj.time(i) <= loadTime )
                      obj.equivForce = obj.force ;
                  end

                  % % 4.2.2. add other parts to equivForce vector. Note initial data are zero
                  if( i > 1 )
                      tempU = ( alpha0 + alpha * alpha1 ) * U + ( alpha2 + alpha * alpha4 ) * Udot + ( alpha3 + alpha * alpha5 ) * Uddot ;
                      obj.equivForce = obj.equivForce + obj.mass * tempU ;

                      tempU =  ( beta0 * alpha1 ) * U + ( beta0 * alpha4 ) * Udot + ( beta0 * alpha5 ) * Uddot ;
                      obj.equivForce = obj.equivForce + obj.stiffness * tempU ;
                  end

                  % % 4.3. Solve equivalent system
                   obj.equivSolution = obj.equivStiffness \ obj.equivForce ;


                  % % 4.4. update velocity and acceleration
                  UdotNew  = alpha1 * ( obj.equivSolution - U ) - alpha4 * Udot - alpha5 * Uddot ;
                  UddotNew = alpha0 * ( obj.equivSolution - U ) - alpha2 * Udot - alpha3 * Uddot ;

                  % % 4.5. Redefine variables
                  U     = obj.equivSolution    ;
                  Udot  = UdotNew              ;
                  Uddot = UddotNew             ;

                  % % 4.6.  displacement vector for the tip function


                  for  n = 1:3
                       var_id    = obj.degreeOfFreedom( obj.tipNodeId  , n ) ;
                       if( var_id == 0 )
                           continue ;
                       end
                       obj.displacement(i+1, n) = U(    var_id ) ;
                  end % of for-n loop

                  x = obj.displacement(i+1, :) ;
                  
                  fprintf('Tip Displacement Vector = (%12.8f , %12.8f , %12.8f ) \n',x(1) , x(2) , x(3) ) ;
                  fprintf('\n\n' ) ;
                           
              end % of for-i loop



              % % 5. Free memory

              U              = [] ;
              Udot           = [] ;
              Uddot          = [] ;
              UdotNew        = [] ;
              UddotNew       = [] ;
              tempU          = [] ;

              obj.mass       = [] ;
              obj.stiffness  = [] ;
              obj.force      = [] ;
              obj.degreeOfFreedom = [] ;
              obj.equivStiffness  = [] ;
              obj.equivForce      = [] ;
              obj.equivSolution   = [] ;



        end % of function doTransientAnalysis


        %%=================================================================================%%
        function  fitToDecayingCurve( obj , geomObj , Load )
            
              % % plot components of displacement vector

              figure()
              subplot(3,1,1);
              plot( obj.time ,obj.displacement(:,1) , 'r' , 'LineWidth' , 2 )
              ylabel(' $u_x$ (meter)' , 'interpreter' , 'latex')
              title('Different components of displacement vector at tip point vs. time')
              set(gca , 'FontSize' , 16)

              subplot(3,1,2);
              plot( obj.time , obj.displacement(:,2) , 'g' , 'LineWidth' , 2 )
              ylabel(' $u_y $ (meter)' , 'interpreter' , 'latex')
              set(gca , 'FontSize' , 16)

              subplot(3,1,3);
              plot( obj.time , obj.displacement(:,3) , 'b' , 'LineWidth' , 2 )
              ylabel(' $u_z $ (meter)' , 'interpreter' , 'latex')
              xlabel('time (sec)')
              set(gca , 'FontSize' , 16)


              %% plot magnitude of displacement vector
              disp_mag = zeros( [size(obj.displacement , 1) ,1] , 'double') ;
              for i = 1 : length( obj.time )
                  disp_mag(i,1) = norm( obj.displacement(i,:) ) ;
              end

              figure()
              plot( obj.time ,disp_mag , 'r' , 'LineWidth' , 2 )
              xlabel('time (sec)')
              ylabel(' $u = \sqrt{ u_x^2 + u_y^2 + u_z^2 }$ (meter)' , 'interpreter' , 'latex')
              title('Magnitude of displacement vector at tip node vs. time')
              set(gca , 'FontSize' , 16)


              %% get envelope
              difMag = diff( disp_mag ) ;
              numMax = 0 ;
              idMax = zeros( size( difMag ) , 'int32')  ;

              for i = 1 : length( difMag ) -1
                  if( ( difMag(i) > 0 ) && ( difMag(i+1) < 0 ) )
                      numMax = numMax + 1 ;
                      idMax(numMax) = i+1 ;
                  end
              end

              if( numMax < 2 )
                  error(' Not enough data for curve fitting: the curve should have at least two maximum.')
              end
              idMax = idMax(1:numMax)     ;

              %% fit envelope to exponential function
              t_env = obj.time( idMax(1:numMax) )    ;
              u_env = disp_mag( idMax(1:numMax) )    ;

              [xData, yData] = prepareCurveData( t_env, u_env );

              % % get initial data.
              zz = log( yData ) ;
              M = [ mean( xData .^2 ) , -mean(xData) ; mean(xData) , -1.0 ] ;
              b = [ -mean( zz .* xData ) ; -mean( zz ) ] ;
              xx = M \b ;

              % Set up fittype and options.
              ft = fittype( 'exp1' );
              opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
              opts.Display = 'Off';
              opts.StartPoint = [exp(xx(2)) , -xx(1)] ;
              opts.TolFun = 1.0e-8 ;
              opts.TolX   = 1.0e-8 ;

              % Fit model to data.
              [fitresult, ~] = fit( xData, yData, ft, opts ) ;

              u_fitted = fitresult.a * exp( fitresult.b * obj.time ) ;

              % Plot fit with data.
              figure( 'Name', 'Magnitude of displacement vs. time and its envelope' );
              plot( obj.time ,disp_mag , 'r' , 'LineWidth' , 2 )
              hold on
              plot( obj.time ,u_fitted , 'b' , 'LineWidth' , 2 )
              hold off
              xlabel('time (sec)')
              ylabel(' $u = \sqrt{ u_x^2 + u_y^2 + u_z^2 }$ (meter)' , 'interpreter' , 'latex')
              legend( 'Magnitude of displacement: $u = \sqrt{ u_x^2 + u_y^2 + u_z^2 }$' , 'Fitted envelope curve: $y(t) = A e^{-\gamma t }$', 'Location', 'NorthEast', 'Interpreter', 'latex' );
              title('Magnitude of displacement vector at tip node vs. time')
              grid on
              set(gca , 'FontSize' , 16)


              % % printing results.
              rtip = geomObj.nodes( obj.tipNodeId , : ) ;
              fileID = fopen('./output/TimeDomainAnalysisInfo.txt','w');
              fprintf(fileID,' Result of Time domain analysis: \n'  );
              fprintf(fileID,'\n'  );
              fprintf(fileID,' Tip position            = \t (%15.12f , %15.12f , %15.12f ) \n' ,  rtip(1) , rtip(2) , rtip(3) );
              fprintf(fileID,' Load acted on Tip       = \t (%15.12f , %15.12f , %15.12f ) \n' ,  Load(1) , Load(2) , Load(3) );
              fprintf(fileID,' Min time step               = \t  %15.12f  \n' ,  obj.time(2)- obj.time(1) );
              fprintf(fileID,' Max time step               = \t  %15.12f  \n' ,  obj.time(length(obj.time))- obj.time(length(obj.time) - 1) );
              fprintf(fileID,' Analysis Time Interval  = \t  [%15.12f , %15.12f]  \n' ,  obj.time(1) , obj.time(end) );

              fprintf(fileID,'\n'  );
              fprintf(fileID,'\n'  );

              fprintf(fileID,' Mesh Information: \n'  );
              fprintf(fileID,' Number of Nodes      = \t %8d     \n', geomObj.numNodes );
              fprintf(fileID,' Number of Elements   = \t %8d     \n', geomObj.numElements );
              fprintf(fileID,' Max Element Size     = \t %15.12d  \n', geomObj.maxElementSize );
              fprintf(fileID,' Min Element Size     = \t %15.12d  \n', geomObj.minElementSize );

              fprintf(fileID,'\n'  );
              fprintf(fileID,'\n'  );
              fprintf(fileID,'\n'  );

              fprintf(fileID,' Max (Magnitude of) Oscillation Amplitude of Beam at tip: \n'  );
              [ max_amplitude , timeStep_id ] = max( disp_mag ) ;
              fprintf(fileID,' At Time Iteration   = \t %15.12f   \n' , obj.time(timeStep_id) );
              fprintf(fileID,' With magnitude      = \t %15.12f   \n' , max_amplitude );


              %% 13. Write fit data

              fprintf(fileID,'\n'  );
              fprintf(fileID,'\n'  );
              fprintf(fileID,'\n'  );

              fprintf(fileID,' Fitting curve to envelope of magnitude of displacement vector to the function: \n'  );
              fprintf(fileID,'       y(t) = A * exp( - gamma * t )  \n'  );
              fprintf(fileID,'\n'  );
              fprintf(fileID,'Fit Result: \n'  );
              fprintf(fileID,'  Ocillation Domain    : A       = \t %15.12f ( meter   ) \n' , fitresult.a );
              fprintf(fileID,'  Decay Factor         : gamma   = \t %15.12f ( 1/sec   ) \n' , fitresult.b );
              fprintf(fileID,"  Decay Rate (Time to reach 37%% of max Amplitude)    = \t %15.12f ( sec )  \n" , 1.0 / abs(fitresult.b)  );
              fprintf(fileID,'\n'  );
              fclose(fileID);

        end % of function fitToDecayingCurve

  end % of public methods


  %% private methods
  methods (Access = private)


        %%=======================================================================================%%
        function doFEMforElement(obj , geomObj , elemId , rho )

            % 1. set elemental matrices initially to zero
            obj.elementalMass(:,:)     = 0.0 ;
            obj.elementalStiffness(:,:) = 0.0 ;

            obj.rowVarSet(:)  = 0  ;
            numRow = 0 ;  % numVars in row var set

            numDoF = 0 ;

            % % 2. set rows and column variable sets
            for i = 1:size(geomObj.elements , 2)
                n = geomObj.elements( elemId , i ) ;
                for j = 1:3  % num dim
                    numRow = numRow + 1 ;
                    obj.rowVarSet(numRow) = obj.degreeOfFreedom(n,j) ;
                end % for
            end % of for i-loop

            % % if the element has no dof then do nothing.
            if( ~any( obj.rowVarSet(1:numRow)) )
                return ;
            end

            % % 3. Calculate elemental mass and stiffness matrices
            % % 3.1. Loop over Gauss points for integrating functions in 3D, i.e. three loops
            for i1 = 1: size( obj.gaussPoint , 1 )
                        % % 3.2. calculate chape function at given point
                        obj.getShapeFunctionAt( obj.gaussPoint(i1,:)   ) ;

                        % % 3.3. calculate local derivative of shape function at given Gauss point
                        obj.getLocalDerShapeFunctionAt( obj.gaussPoint(i1,:)  ) ;

                        % % 3.3. calculate local derivative of shape function at given Gauss point
                        obj.getJacobainMatrix( geomObj , elemId ) ;

                        % % 3.4. calculate global derivative of shape function at given Gauss point
                        obj.getGlobalDerShapeFunction() ;

                        % % 3.5. calculate volume element
                        detJacobian = abs( det( obj.JacobianMatrix ) ) ;

                        if( detJacobian < 1.0e-14)
                            error( strcat( 'Jacobian matrix is singular for element = ' , num2str(elemId)  ,...
                                           ', at Gauss pointid  =' , num2str(i1)  )) ;
                        end

                        dVolume = obj.gaussWeight(i1) * detJacobian  ;

                        % % 3.6. calculate elemental mass and stiffness matrices
                        for row = 1 : 10

                            % form B_{\mu} matrix
                            obj.BmatrixMu(:,:) = 0.0 ;
                            obj.BmatrixMu(1,1) = obj.dGlobalShapeFuncs(1,row) ; % D_x N_{\mu}
                            obj.BmatrixMu(2,2) = obj.dGlobalShapeFuncs(2,row) ; % D_y N_{\mu}
                            obj.BmatrixMu(3,3) = obj.dGlobalShapeFuncs(3,row) ; % D_z N_{\mu}
                            obj.BmatrixMu(4,1) = obj.BmatrixMu(2,2) ;           % D_y N_{\mu}
                            obj.BmatrixMu(4,2) = obj.BmatrixMu(1,1) ;           % D_x N_{\mu}
                            obj.BmatrixMu(5,2) = obj.BmatrixMu(3,3) ;           % D_z N_{\mu}
                            obj.BmatrixMu(5,3) = obj.BmatrixMu(2,2) ;           % D_y N_{\mu}
                            obj.BmatrixMu(6,1) = obj.BmatrixMu(3,3) ;           % D_z N_{\mu}
                            obj.BmatrixMu(6,3) = obj.BmatrixMu(1,1) ;           % D_x N_{\mu}


                            for col = 1:10

                                mVal = rho * obj.shapeFunctionVector(row) * obj.shapeFunctionVector( col ) * dVolume ;
                                r1 = 3 * ( row - 1 ) ;  % begin of row block
                                c1 = 3 * ( col - 1 ) ;  % begin of col block

                                % % 3.6.2. add to mass matrix, it is a multiple of I_3, the identity matrix
                                for n = 1 : 3  % dimension of structure
                                    obj.elementalMass(r1 + n ,c1 + n) = obj.elementalMass(r1 + n ,c1 + n) + mVal ;
                                end % for n-loop

                                % % 3.6.3. Form B matrix
                                % form B_{\mu} matrix
                                obj.BmatrixNu(:,:) = 0.0 ;
                                obj.BmatrixNu(1,1) = obj.dGlobalShapeFuncs(1,col) ; % D_x N_{\nu}
                                obj.BmatrixNu(2,2) = obj.dGlobalShapeFuncs(2,col) ; % D_y N_{\nu}
                                obj.BmatrixNu(3,3) = obj.dGlobalShapeFuncs(3,col) ; % D_z N_{\nu}
                                obj.BmatrixNu(4,1) = obj.BmatrixNu(2,2) ;           % D_y N_{\nu}
                                obj.BmatrixNu(4,2) = obj.BmatrixNu(1,1) ;           % D_x N_{\nu}
                                obj.BmatrixNu(5,2) = obj.BmatrixNu(3,3) ;           % D_z N_{\nu}
                                obj.BmatrixNu(5,3) = obj.BmatrixNu(2,2) ;           % D_y N_{\nu}
                                obj.BmatrixNu(6,1) = obj.BmatrixNu(3,3) ;           % D_z N_{\nu}
                                obj.BmatrixNu(6,3) = obj.BmatrixNu(1,1) ;           % D_x N_{\nu}


                                kblock = transpose( obj.BmatrixMu ) * obj.rigidity * obj.BmatrixNu ;
                                kblock = kblock * dVolume ;

                                obj.elementalStiffness(r1+1:r1+3 , c1+1:c1+3) = obj.elementalStiffness(r1+1:r1+3 , c1+1:c1+3) + kblock ;

                            end % of for-col loop
                        end  % of for-row  loop
            end % of for-i1 loop

            % % 4. Add to Global matrices
            obj.addToMassMatrix( numRow ) ;
            obj.addToStiffnessMatrix( numRow ) ;


        end % of function doFEMforElement

        %%=======================================================================================%%
        function getShapeFunctionAt(obj, a_point  )
            xi   = a_point(1) ;
            eta  = a_point(2) ;
            zeta = a_point(3)  ;
            t = [ 1.0 - xi - eta - zeta , xi , eta , zeta ] ;
            obj.shapeFunctionVector = [ t(1) * ( 2.0 * t(1) - 1.0 ) , t(2) * ( 2.0 * t(2) - 1.0 ) ,...
                                        t(3) * ( 2.0 * t(3) - 1.0 ) , t(4) * ( 2.0 * t(4) - 1.0 ) ,...
                                        4.0 * t(1) * t(2) , 4.0 * t(2) * t(3) , 4.0 * t(1) * t(3) ,...
                                        4.0 * t(1) * t(4) , 4.0 * t(2) * t(4) , 4.0 * t(3) * t(4) ] ;
        end % of function getShapeFunction


        %%=======================================================================================%%
        function getLocalDerShapeFunctionAt(obj , a_point )

             xi   = a_point(1) ;
             eta  = a_point(2) ;
             zeta = a_point(3) ;

             obj.dLocalShapeFuns(:,:) = 0.0 ;
             obj.dLocalShapeFuns( 1:3 , 1 ) =  4.0 *( xi  + eta + zeta ) - 3.0 ;
             obj.dLocalShapeFuns( 1   , 2 ) =  4.0 * xi   - 1.0 ;
             obj.dLocalShapeFuns( 2   , 3 ) =  4.0 * eta  - 1.0 ;
             obj.dLocalShapeFuns( 3   , 4 ) =  4.0 * zeta - 1.0 ;

             obj.dLocalShapeFuns( 1   , 5 ) =  4.0 *( 1.0 - 2.0 *xi  - eta - zeta )  ;
             obj.dLocalShapeFuns( 2:3 , 5 ) = -4.0 * xi  ;

             obj.dLocalShapeFuns( 1   , 6 ) =  4.0 * eta  ;
             obj.dLocalShapeFuns( 2   , 6 ) =  4.0 * xi   ;

             obj.dLocalShapeFuns( 1   , 7 ) = -4.0 * eta  ;
             obj.dLocalShapeFuns( 2   , 7 ) =  4.0 *( 1.0 - xi  - 2.0 * eta - zeta )  ;
             obj.dLocalShapeFuns( 3   , 7 ) = -4.0 * eta  ;

             obj.dLocalShapeFuns( 1:2 , 8 ) = -4.0 * zeta ;
             obj.dLocalShapeFuns( 3   , 8 ) =  4.0 * ( 1.0 - xi  -  eta - 2.0 * zeta )  ;

             obj.dLocalShapeFuns( 1   , 9 ) =  4.0 * zeta ;
             obj.dLocalShapeFuns( 3   , 9 ) =  4.0 * xi   ;

             obj.dLocalShapeFuns( 2   , 10) =  4.0 * zeta ;
             obj.dLocalShapeFuns( 3   , 10) =  4.0 * eta  ;

        end % of function getLocalDerShapeFunctionAt

        %%=======================================================================================%%
        function  getJacobainMatrix( obj , geomObj, elemId )
            obj.JacobianMatrix(:,:) = 0.0 ;

            for i = 1 : 3
                for j = 1 : 3
                    for n = 1 :10   % % loop over nodes
                        obj.JacobianMatrix(i,j) = obj.JacobianMatrix(i,j) + obj.dLocalShapeFuns(j,n) *...
                                                  geomObj.nodes( geomObj.elements( elemId , n) , i ) ;
                    end
                end
            end

        end % of function  getJacobainMatrixAt

        %%=====================================================================================%%
        function  getGlobalDerShapeFunction( obj )

            obj.inverseJacobianT  = inv( transpose( obj.JacobianMatrix ) ) ;
            obj.dGlobalShapeFuncs = obj.inverseJacobianT * obj.dLocalShapeFuns ;

        end % of function getGlobalDerShapeFunctionAt

        %%=====================================================================================%%
        function addToMassMatrix( obj , numRow )


             for i = 1 : numRow
                 if( obj.rowVarSet(i) == 0 )
                     continue ;
                 end

                 for j = 1 : numRow
                     if( obj.rowVarSet(j) == 0 )
                         continue ;
                     end
                     obj.mass( obj.rowVarSet(i) ,obj.rowVarSet(j) ) = obj.mass( obj.rowVarSet(i) ,obj.rowVarSet(j) ) + obj.elementalMass(i,j) ;

                 end % of for-j loop
             end % of for-i loop

        end % of function addToMassMatrix

        %%=====================================================================================%%
        function addToStiffnessMatrix( obj , numRow )

           for i = 1 : numRow
               if( obj.rowVarSet(i) == 0 )
                   continue ;
               end

               for j = 1 : numRow
                   if( obj.rowVarSet(j) == 0 )
                       continue ;
                   end
                   obj.stiffness( obj.rowVarSet(i) ,obj.rowVarSet(j) ) = obj.stiffness( obj.rowVarSet(i) ,obj.rowVarSet(j) ) + obj.elementalStiffness(i,j) ;

               end % of for-j loop
           end % of for-i loop
        end % of function addToStiffnessMatrix


  end % of private methods

end % of class




