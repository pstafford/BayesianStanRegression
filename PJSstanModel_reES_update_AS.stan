functions {
    // define the functional form based upon Abrahamson & Silva (2014)
    // note that this differs from PJSstanModel_reES_AS.stan in that here a fixed effects vector is passed
    real PJSgmmAS14(real m, real r, real v, real Fnm, real Frv, vector bi, real eta, real deltaS2S) {
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
        real a1 = bi[1] + eta + deltaS2S; // put the random effects for event & site here on the constant coefficient so it propagates through site response
        real a2 = bi[2];
        real a3 = bi[3];
        real a4 = bi[4];
        real a5 = bi[5];
        real a6 = bi[6];
        real a8 = bi[7];
        real a10 = bi[8];
        real a11 = bi[9];
        real a12 = bi[10];
        real a17 = bi[11];
        real b = bi[12];

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
            f1 = a1 + a5*(m-M1) + a8*(8.5-m)^2 + (a2 + a3*(m-M1))*log(R) + a17*r;
        } else if ( m < M2 ) {
            f1 = a1 + a4*(M2-M1) + a8*(8.5-M2)^2 + a6*(m-M2) + (a2 + a3*(M2-M1))*log(R) + a17*r;
        } else {
            f1 = a1 + a4*(m-M1) + a8*(8.5-m)^2 + (a2 + a3*(m-M1))*log(R) + a17*r;
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
    int<lower=0> N;                     // number of new observations
    int<lower=1> I;                     // total number of events (previous & current)
    int<lower=1> Ip;                    // number of events already used in prior
    int<lower=1> J;                     // number of stations (including future stations)
    int<lower=1> K;                     // number of fixed effects
    vector[N] y;                        // logarithmic intensity measures
    real M;                             // new magnitude value
    real Rrup[N];                       // new distance values
    real Vs30[N];                       // new shear-wave velocities
    real Fnm;                           // new normal mechanism flag
    real Frv;                           // new reverse mechanism flag
    int<lower=1,upper=I> EventID;       // new event identifier
    int<lower=1,upper=J> StationID[N];  // station identifier
    int<lower=0,upper=1> NewStation[J]; // is this a new station (1 = yes, 0 = no)
    vector[K] MuFixed;                  // prior estimates of the fixed effects
    vector[K] SigmaFixed;               // error estimates for the fixed effects
    corr_matrix[K] OmegaFixed;          // correlation for prior on fixed effects
    real MuTau;                         // current tau estimate
    real SigmaTau;                      // current error in tau
    real MuPhiS2S;                      // current phiS2S estimate
    real SigmaPhiS2S;                   // current error in phiS2S
    real MuPhi;                         // current phi estimate
    real SigmaPhi;                      // current error in phi
    real MuRandomEta[Ip];               // current estimates of random effects (for priors)
    real SigmaRandomEta[Ip];            // error estimates in the random effects
    real MuRandomDelta[J];              // current estimates of random effects (for priors)
    real SigmaRandomDelta[J];           // error estimates in the random effects
}
transformed data {
    cov_matrix[K] CovFixed = quad_form_diag( OmegaFixed, SigmaFixed );
}
parameters {
    vector[K] bi;                       // fixed effects regression parameters
    real<lower=0> tau;                  // inter-event variability
    real<lower=0> phiS2S;               // inter-site variability
    real<lower=0> phi;                  // intra-event variability
    vector[I] eta;                      // inter-event residuals (all events)
    vector[J] delta;                    // inter-station residuals (all stations)
}
model {
    // define strongly informative priors from the previous model run
    bi ~ multi_normal( MuFixed, CovFixed );
    tau ~ normal( MuTau, SigmaTau );
    phiS2S ~ normal( MuPhiS2S, SigmaPhiS2S );
    phi ~ normal( MuPhi, SigmaPhi );

    // note that the last entries of MuRandom and SigmaRandom should be 0, tau
    // other entries are the previous estimates and errors for the random effects
    // Note that we update the past event terms here to reflect changes in the fixed effects
    // if the new model intercept (bi[1]) is greater than the previous (MuFixed[1]), then
    // we need to reduce the magnitudes of the historical random effects
    for ( i in 1:I ) {
        if ( i < I ) {
            eta[i] ~ normal( MuRandomEta[i] - (bi[1] - MuFixed[1]), SigmaRandomEta[i] );
        } else {
            eta[i] ~ normal( 0, tau );
        }
    }
    // need more care here to adjust only previously seen stations
    for ( i in 1:J ) {
        if ( NewStation[i] == 1 ) {
            delta[i] ~ normal( 0, phiS2S );
        } else {
            delta[i] ~ normal( MuRandomDelta[i] - (bi[1] - MuFixed[1]), SigmaRandomDelta[i] );
        }
    }

    for ( i in 1:N ) {
        y[i] ~ normal( PJSgmmAS14(M, Rrup[i], Vs30[i], Fnm, Frv, bi, eta[EventID], delta[StationID[i]]), phi );
    }

}
