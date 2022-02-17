function c_ts = calculateCorrelationFunction(weights,omegas,ts,beta,Delta_E)

    c_ts = 1.0i* Delta_E*ts-sum(((weights./omegas)').* (coth(0.5*beta*omegas').*(1-cos(omegas'.*ts))+1.0i * sin(omegas'.*ts)),1) ;
    c_ts = exp(c_ts) ;

end