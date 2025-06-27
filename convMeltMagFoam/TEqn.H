{
    int iter = 0;
    scalar residual = 1;
    scalar meanResidual = 1;

    do
    {
        iter++;
        alpha1.storePrevIter();

        cp = alpha1*cpl + alpha2*cps;
        phiCp = fvc::interpolate(cp)*phi;
        lambda = alpha1*lambdal + alpha2*lambdas;

        fvScalarMatrix TEqn
        (
            fvm::ddt(cp, T)
          + fvm::div(phiCp, T)
          + L*fvc::ddt(alpha1)
          + L*fvc::div(phi, alpha1)
          - Tmelt*fvc::ddt(cp)
          - Tmelt*fvc::div(phiCp)
          - fvm::laplacian(lambda/rho, T)
        );

        TEqn.solve();

        Tcorr = (Tl-Ts)*alpha1+Ts;  
        alpha1 = max(min(alpha1+alphaRel*cp/L*(T-Tcorr),scalar(1)),scalar(0));
        alpha2 = scalar(1) - alpha1;
        
        residual = max(mag(alpha1.internalField()-alpha1.prevIter().internalField())).value();
        meanResidual = (sum(mag(alpha1.internalField() - alpha1.prevIter().internalField()) * mesh.V()).value()) / sum(mesh.V()).value();
   
        //meanResidual = sum(mag(alpha1.internalField()-alpha1.prevIter().internalField())*mesh.V())/sum(mesh.V()).value(); 


        /*
        residual = max(mag(alpha1.internalField()-alpha1.prevIter().internalField()));
        meanResidual = sum(mag(alpha1.internalField()-alpha1.prevIter().internalField())*mesh.V())/sum(mesh.V()).value();
          */
        
        Info << "Correcting alpha1, mean residual = " << meanResidual
             << ", max residual = " << residual
             << endl;

    } while ((iter < minTCorr || residual > alphaTol ) && iter <= maxTCorr);

    h = cp*(T-Tmelt)+alpha1*L;
    DC = DCl*pow(alpha2,2)/(pow(alpha1,3)+DCs);
    rhok = scalar(1)-max(beta*(T-Tl),scalar(0));
}
