







tmp<volScalarField> wetSteamSystem::L() const
{

    volScalarField Ts = TSat();

    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        {return gasProps_->rho(pp,TT);};
        
    volScalarField rhos = applyFunction2(f, gas().p(), Ts, "rhos", dimDensity);

    return tmp<volScalarField>
        (
	    //new volScalarField("L", saturation_->dpsdT(Ts)*Ts/rhos)
            new volScalarField("L", gas().Cp()*(T() - Ts) + saturation_->dpsdT(Ts)*Ts/rhos)
        );

    /*const std::function<scalar(scalar)> f = [&](scalar TT)
        {return 461.52*(-2.7246e-2*sqr(TT) + 2*1.6853e-5*pow(TT,3) + 2.4576*TT + 6094.4642);};

    return applyFunction1(f, T(), "L", dimEnergy/dimMass);*/

}