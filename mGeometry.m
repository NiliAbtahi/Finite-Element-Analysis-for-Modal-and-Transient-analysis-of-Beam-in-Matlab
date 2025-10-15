%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The mGeometry is a handle class (in contrast to a value class) to avoid
% unnecessary copying. It is the matlab way for working with class
% reference.
%
% It contains a standard pdeModel object from pde toolbox. This object is used to
% access reading STL file, rescaling geometry and meshing it.
% No more application of pde model will be used.
% Once the mesh is generated the pdeModel will be deleted to free memory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef  mGeometry < handle

    % % public properties==========================================================================

    properties
       pdeModel        = []  ;  % pde model: required only for reading, rescaling and meshing geometry
       nodes           = []  ;  % mesh nodes
       elements        = []  ;  % mesh elements
       maxElementSize  = 0   ; % max elem length in the mesh
       minElementSize  = 0   ; % min elem length in the mesh
       tip_position    = []  ; % position of tip vertex where the load is axerted.
       numNodes        = 0   ; % number of nodes in the mesh
       numElements     = 0   ; % number of elements in the mesh
       fix_face_axis   = 0   ; % The fixed face axis is one of the items: 1 (xmin or xmax), 2(ymin or ymax), 3(zmin or zmax)
       fix_face_value  = 0.0 ; % The value specifying fixed plane in 3D like xmin = val, xmax = val, etc.

    end  % of public properties


    % % public methods
    methods

        %%=======================================================================================%%
        function obj = mGeometry( stl_file_name ) %% constructor
            if nargin == 1
               obj.pdeModel    = createpde;
               obj.pdeModel.Geometry = importGeometry(obj.pdeModel , stl_file_name) ;
            else
                error("The constructor of geometry class needs address of stl file as input.")
            end
        end


        %%=======================================================================================%%
        function  rescaleGeometry(obj , scaleVector )

            % 1. Get vertices of the geometry and find Length, width and height of current geometry
            V = obj.pdeModel.Geometry.vertexCoordinates(1:obj.pdeModel.Geometry.NumVertices) ;
            xmin = min( V(:,1) ) ;
            xmax = max( V(:,1) ) ;
            ymin = min( V(:,2) ) ;
            ymax = max( V(:,2) ) ;
            zmin = min( V(:,3) ) ;
            zmax = max( V(:,3) ) ;

            curr_length = xmax - xmin ;
            curr_width  = ymax - ymin ;
            curr_height = zmax - zmin ;

            % % 1.1.2. Get vector of scale factor, i.e. [s_x ,s_y , s_z ] where:
            % %        s_x is the ratio of user-defined Length and that of stl file
            % %        s_y is the ratio of user-defined Width  and that of stl file
            % %        s_z is the ratio of user-defined Heigth and that of stl file

            rescale_vector = [ scaleVector(1) / curr_length , scaleVector(2) / curr_width , scaleVector(3) / curr_height ] ;

            obj.pdeModel.Geometry = scale(obj.pdeModel.Geometry , rescale_vector ) ;



        end

        %%======================================================================================%%
        function  specifyTipVertex(obj)

            % 1. Get vertices of the geometry and find Length, width and height of current geometry
            % % 1.1. specify tip position, i.e. the point at which the force ( or impulse/load ) is
            % %      exerted.

            % % 1.2. Get vertices of the rescaled geometry.

            V = obj.pdeModel.Geometry.vertexCoordinates(1:obj.pdeModel.Geometry.NumVertices) ;

            % % 1.3. This point is selected at the free end of the beam, i.e. at the plane x = xmin,


            xmin = min( V(:,1) )  ;
            ymin = min( V(:,2) )  ;
            ymax = max( V(:,2) )  ;

            % % 1.4. find the point at line x = xmin , y = ( ymin + ymax ) /2 and max(z) ;
            mask_x = abs( V(:,1) - xmin ) < 1.0e-3 ;
            [ zmax ,~  ] = max( V(mask_x,3) ) ;
            obj.tip_position = [ xmin  , (ymin + ymax)/ 2.0 , zmax ] ;

            figure()
            pdegplot(obj.pdeModel)
            hold on
               plot3( obj.tip_position(1) , obj.tip_position(2) , obj.tip_position(3) ,...
                     'o','Color','r','MarkerSize',15,'MarkerFaceColor','r') ;
            hold off
            title('Rescaled Geometry and Tip Vertex')
            set(gca , 'FontSize' , 18)


        end


        %%=======================================================================================%%
        function createMesh(obj)
            generateMesh( obj.pdeModel );

            obj.nodes          = transpose( obj.pdeModel.Mesh.Nodes    ) ;
            obj.elements       = transpose( obj.pdeModel.Mesh.Elements ) ;
            obj.maxElementSize = obj.pdeModel.Mesh.MaxElementSize ;
            obj.minElementSize = obj.pdeModel.Mesh.MinElementSize ;


            figure()
            pdeplot3D( obj.pdeModel )
            title('Mesh')
            set(gca , 'FontSize' , 18)

            % % free arrays
            obj.pdeModel = [] ;

            obj.numNodes    = size( obj.nodes    , 1 );
            obj.numElements = size( obj.elements , 1 );
        end

        %%====================================================================================%%
        function fixFaceAt( obj , plane_location )

            if( strcmp( plane_location , 'xmax' ) )
                obj.fix_face_axis = 1 ;
                obj.fix_face_value = max( obj.nodes(:,1) ) ;

            elseif( strcmp( tip_surface_location , 'xmin' ) )
                obj.fix_face_axis = 1 ;
                obj.fix_face_value = min( obj.nodes(:,1) ) ;

             elseif( strcmp( tip_surface_location , 'ymax' ) )
                obj.fix_face_axis = 2 ;
                obj.fix_face_value = max( obj.nodes(:,2) ) ;

             elseif( strcmp( tip_surface_location , 'ymin' ) )
                obj.fix_face_axis = 2 ;
                obj.fix_face_value = min( obj.nodes(:,2) ) ;

             elseif( strcmp( tip_surface_location , 'zmax' ) )
                obj.fix_face_axis = 3 ;
                obj.fix_face_value = max( obj.nodes(:,3) ) ;

             elseif( strcmp( tip_surface_location , 'zmin' ) )
                obj.fix_face_axis = 3 ;
                obj.fix_face_value = min( obj.nodes(:,3) ) ;
             else
                error("The argument of fixFace function should be one of: xmin,xmax,ymin,ymax,zmin and zmax")
            end
        end

        %%====================================================================================%%
        function writeGeometryData( obj )

           % 1. write nodes to node.txt file
           fileID = fopen('./meshData/node.txt','w');
           for i = 1: obj.numNodes
               fprintf(fileID,'%12.8f ,%12.8f ,%12.8f \n',obj.nodes(i,:) );
           end
           fclose(fileID);

           % 2. write elements to element.txt file
           fileID = fopen('./meshData/element.txt','w');
           for i = 1: obj.numElements
               fprintf(fileID,'%d , %d , %d , %d , %d , %d , %d , %d , %d , %d \n',obj.elements(i,:) );
           end
           fclose(fileID);


            % 3. write other data to other.txt file
           fileID = fopen('./meshData/other.txt','w');
           fprintf(fileID,'%d , %d \n',obj.numNodes , obj.numElements );
           fprintf(fileID,'%12.8f ,%12.8f ,%12.8f \n',obj.tip_position(1:3) );
           fclose(fileID);
        end


        %%====================================================================================%%
        function  destrcutGeometry(obj)
            obj.pdeModel        = []  ;  % pde model: required only for reading, rescaling and meshing geometry
            obj.nodes           = []  ;  % mesh nodes
            obj.elements        = []  ;  % mesh elements
            obj.tip_position    = []  ; % position of tip vertex where the load is axerted.
            obj = [] ;
        end % of destructor function


    end % of public methods

end % of class
