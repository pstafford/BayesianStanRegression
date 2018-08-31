functions {
    // define the functional form based upon Abrahamson & Silva (2014)
    real PJSgmmAS14(real m, real r, real v, int Fnm, int Frv, real a1, real a2, real a3, real a4, real a5, real a6, real a8, real a10, real a11, real a12, real a17, real b, real eta, real deltaS2S) {
        // declare some variables
        real R;
        real f1;
        real f5;
        real lnSa;

        // parameters constrained externally
        real V1 = 1500.0; // this is assuming T <= 0.5 seconds
        real c = 2.4;
        real n = 1.5;
        real Vlin = 660.0;
        real Vs30star = ( v < V1 ) ? v : V1;

        // extract the coefficients
        real a1e = a1 + eta + deltaS2S; // put the random effect for event and linear site response here on the constant coefficient so it propagates through site response

        // parameters constrained externally
        real M1 = 6.75;
        real M2 = 5.0;
        real c4 = 4.5;

        // declare parameters
        real f7 = 0.0;
        real f8 = 0.0;

        // fictitious depth
        real c4m;
        if ( m > 5.0 ){
            c4m = c4;
        } else if ( m <= 4.0 ) {
            c4m = 1.0;
        } else {
            c4m = c4 - (c4-1.0)*(5.0-m);
        }

        // general distance
        R = sqrt( r*r + c4m*c4m );
        // basic magnitude and distance scaling
        if ( m > M1 ) {
            f1 = a1e + a5*(m-M1) + a8*(8.5-m)^2 + (a2 + a3*(m-M1))*log(R) + a17*r;
        } else if ( m < M2 ) {
            f1 = a1e + a4*(M2-M1) + a8*(8.5-M2)^2 + a6*(m-M2) + (a2 + a3*(M2-M1))*log(R) + a17*r;
        } else {
            f1 = a1e + a4*(m-M1) + a8*(8.5-m)^2 + (a2 + a3*(m-M1))*log(R) + a17*r;
        }
        // style of faulting
        if ( m > 5.0 ) {
            f7 = a11;
            f8 = a12;
        } else if ( m >= 4.0 ) {
            f7 = a11*(m-4.0);
            f8 = a12*(m-4.0);
        }
        // site response
        if ( v >= Vlin ) {
            // linear Vs30 dependent response (independent above V1)
            f5 = (a10 + b*n)*log(Vs30star/Vlin);
        } else {
            // nonlinear response
            // compute Sa1180
            real vnl = 1180.0;
            real lnSa1180 = f1 + Fnm*f7 + Frv*f8 + (a10 + b*n)*log((vnl < V1 ? vnl : V1)/Vlin);
            real Sa1180 = exp(lnSa1180);
            f5 = a10*log(Vs30star/Vlin) - b*log(Sa1180 + c) + b*log(Sa1180 + c*(Vs30star/Vlin)^n );
        }
        // combined motion
        lnSa = f1 + Fnm*f7 + Frv*f8 + f5;
        return lnSa;
    }
}
data {
    int<lower=0> N;                     // number of observations
    int<lower=1> I;                     // number of events
    int<lower=1> J;                     // number of stations
    vector[N] y;                        // logarithmic intensity measures
    real M[I];                          // magnitude values
    real Rrup[N];                       // distance values
    real Vs30[J];                       // shear-wave velocities
    int Fnm[I];                         // normal mechanism flag
    int Frv[I];                         // reverse mechanism flag
    int<lower=1,upper=I> EventID[N];    // event identifier
    int<lower=1,upper=J> StationID[N];  // station identifier
}
parameters {
    // fixed effects
    real a1;
    real a2;
    real a3;
    real a4;
    real a5;
    real a6;
    real a8;
    real a10;
    real a11;
    real a12;
    real<upper=0> a17;
    real b;
    real<lower=0> tau;                  // inter-event variability
    real<lower=0> phiS2S;               // inter-site variability
    real<lower=0> phi;                  // intra-event variability
    vector[I] eta;                      // inter-event residuals (all events)
    vector[J] deltaS2S;                 // inter-site residuals (all stations)
}
model {
    // define informative priors
    a1 ~ cauchy( 1.0, 1.0 );
    a2 ~ cauchy( -0.79, 0.5 );
    a3 ~ cauchy( 0.275, 0.5 );
    a4 ~ cauchy( -0.1, 0.3 );
    a5 ~ cauchy( -0.41, 0.5 );
    a6 ~ cauchy( 2.154, 0.5 );
    a8 ~ cauchy( -0.015, 0.3 );
    a10 ~ cauchy( 1.735, 0.5 );
    a11 ~ cauchy( 0.0, 0.3 );
    a12 ~ cauchy( -0.1, 0.3);
    a17 ~ cauchy( -0.0072, 0.05 );
    b ~ cauchy( -1.47, 0.5 );
    tau ~ cauchy( 0.3, 0.5 );
    phiS2S ~ cauchy( 0.3, 0.5 );
    phi ~ cauchy( 0.45, 0.5 );

    eta ~ normal( 0, tau );
    deltaS2S ~ normal( 0, phiS2S );

    for ( i in 1:N ) {
        y[i] ~ normal( PJSgmmAS14(M[EventID[i]], Rrup[i], Vs30[StationID[i]], Fnm[EventID[i]], Frv[EventID[i]], a1, a2, a3, a4, a5, a6, a8, a10, a11, a12, a17, b, eta[EventID[i]], deltaS2S[StationID[i]]), phi );
    }
}
