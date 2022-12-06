        
        //list with all the slides x coordinates 0.0001 0.001
        List<scalar> x_slice = {0.0001, 0.002, 0.004, 0.006};

        forAll(x_slice,i)
        {
        
        //Y LINES
        //-------------------------------------------------------------------------
        scalar pointsinline = 2000;  //number of points in the line at a specific y coordinate
        scalar x_coordinate_slide = x_slice[i]; 
        List<scalar> resultsiny(pointsinline);
	    scalar min_y = mesh.bounds().min().y()/2.0;
        scalar max_y = mesh.bounds().max().y()/2.0;

        
        for (int j =0; j< pointsinline; j++)
        {
        scalar min_y = mesh.bounds().min().y()/2.0;
        scalar max_y = mesh.bounds().max().y()/2.0;
        scalar ycoordinate = min_y + (max_y - min_y)*j/(pointsinline-1);
        
	    point linepoint1(x_coordinate_slide, ycoordinate, 0.0);
        point linepoint2(x_coordinate_slide, ycoordinate, 0.0);
        linepoint2.z() = mesh.bounds().max().z();
        linepoint1.z() = mesh.bounds().min().z();

        int numpoints = 2000;
        List<point> points(numpoints);
        points[0] = linepoint1;
        points[numpoints-1] = linepoint2;
        for (int i = 1; i < numpoints-1; i++)
        {
            points[i].x() = x_coordinate_slide;
            points[i].y() = ycoordinate;
            points[i].z() = mesh.bounds().min().z() + (mesh.bounds().max().z()-mesh.bounds().min().z())*i/(numpoints-1);
         
        }

        List<scalar> rhoPoints(numpoints);
        forAll(points,i)
        {
            label celli = mesh.findCell(points[i]);
            if (celli != -1)
            {
                rhoPoints[i] = rho[celli]*Y[celli];
             //   Info << rhoPoints[i] << " " << "z" << points[i].z() << " " << "y" <<  points[i].y() <<  endl;
            }else
            {
                rhoPoints[i] = 0.0;
            }
        }
	//getchar();

        scalar rhoSum = 0.0;
         for (int i = 0; i < rhoPoints.size() - 1; i++)
        {
            rhoSum += (rhoPoints[i] + rhoPoints[i + 1]) / 2.0 * (points[i + 1].z() - points[i].z()); //trapezoidal rule
        }
        resultsiny[j] = rhoSum;
       // Info << rhoSum << endl;
        }
        
         std::ofstream myfile;
         //use name of slice as file name
        if (std::ifstream( "y_slice_" + std::to_string(x_coordinate_slide) + ".dat"))
        {
            myfile.open("y_slice_" + std::to_string(x_coordinate_slide) + ".dat", std::ios::app);
            myfile << "ZONE T = \"rho\", I = " << pointsinline << ", F = POINT" << std::endl;
            for (int i = 0; i < pointsinline; i++)
            {
               myfile <<  min_y + (max_y - min_y)*i/(pointsinline-1) << " " << resultsiny[i] << std::endl;
            }
        }
        else
        {
            myfile.open("y_slice_" + std::to_string(x_coordinate_slide) + ".dat");
            myfile << "TITLE = \"rho y\"" << std::endl;
            myfile << "VARIABLES = \"y\", \"rho\"" << std::endl;
            myfile << "ZONE T = \"rho\", I = " << pointsinline << ", F = POINT" << std::endl;
            for (int i = 0; i < pointsinline; i++)
            {
                myfile <<  min_y + (max_y - min_y)*i/(pointsinline-1)  << " " << resultsiny[i] << std::endl;
            }
        }
        myfile.close();
        //-------------------------------------------------------------------------

        //Z LINES
        //-------------------------------------------------------------------------
        scalar pointsinlinez = 2000;  //number of points in the line at a specific y coordinate
        scalar x_coordinate_slidez = x_slice[i];
        List<scalar> resultsinz(pointsinlinez);
        scalar min_z = mesh.bounds().min().z()/2.0;
        scalar max_z = mesh.bounds().max().z()/2.0;

        for (int j =0; j< pointsinlinez; j++)
        {
        scalar min_z = mesh.bounds().min().z()/2.0;
        scalar max_z = mesh.bounds().max().z()/2.0;
        scalar zcoordinate = min_z + (max_z - min_z)*j/(pointsinlinez-1);

        point linepoint1(x_coordinate_slidez, 0.0, zcoordinate);
        point linepoint2(x_coordinate_slidez, 0.0, zcoordinate);
        linepoint2.y() = mesh.bounds().max().y();
        linepoint1.y() = mesh.bounds().min().y();

        int numpoints = 2000;
        List<point> points(numpoints);
        points[0] = linepoint1;
        points[numpoints-1] = linepoint2;
        for (int i = 1; i < numpoints-1; i++)
        {
            points[i].x() = x_coordinate_slidez;
            points[i].z() = zcoordinate;
            points[i].y() = mesh.bounds().min().y() + (mesh.bounds().max().y()-mesh.bounds().min().y())*i/(numpoints-1);
        }

        List<scalar> rhoPoints(numpoints);
        forAll(points,i)
        {
            label celli = mesh.findCell(points[i]);
            if (celli != -1)
            {
                rhoPoints[i] = rho[celli]*Y[celli];
            }else
            {
                rhoPoints[i] = 0.0;
            }
        }
        //getchar();

        scalar rhoSum = 0.0;
         for (int i = 0; i < rhoPoints.size() - 1; i++)
        {
            rhoSum += (rhoPoints[i] + rhoPoints[i + 1]) / 2.0 * (points[i + 1].y() - points[i].y()); //trapezoidal rule
        }
        resultsinz[j] = rhoSum;
        //Info << rhoSum << endl;
        }
        
         std::ofstream myfilez;
         //use name of slice as file name
        if (std::ifstream( "z_slice_" + std::to_string(x_coordinate_slidez) + ".dat"))
        {
            myfilez.open("z_slice_" + std::to_string(x_coordinate_slidez) + ".dat", std::ios::app);
            myfile << "ZONE T = \"rho\", I = " << pointsinlinez << ", F = POINT" << std::endl;
            for (int i = 0; i < pointsinlinez; i++)
            {
               myfilez <<  min_z + (max_z - min_z)*i/(pointsinlinez-1) << " " << resultsinz[i] << std::endl;
            }
        }
        else
        {
            myfilez.open("z_slice_" + std::to_string(x_coordinate_slidez) + ".dat");
            myfilez << "TITLE = \"rho z\"" << std::endl;
            myfilez << "VARIABLES = \"z\", \"rho\"" << std::endl;
            for (int i = 0; i < pointsinlinez; i++)
            {
                myfilez <<  min_z + (max_z - min_z)*i/(pointsinlinez-1)  << " " << resultsinz[i] << std::endl;
            }
        }
        myfilez.close();
}
        //-------------------------------------------------------------------------
  return 0;
