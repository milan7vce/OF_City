 {
     scalar sumLocalContErr = runTime.deltaTValue()*

     mag(rhokRK)().weightedAverage(mesh.V()).value();
 
     scalar globalContErr = runTime.deltaTValue()*

     rhokRK.weightedAverage(mesh.V()).value();

     cumulativeContErr += globalContErr;
 
     resFile << runTime.timeName() << "    " << sumLocalContErr << "    "  <<  globalContErr  << "    "  << cumulativeContErr << endl;
 }

