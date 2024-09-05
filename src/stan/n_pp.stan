
    data{
      int k;
      vector[k] y0;
      vector[k] sigma;
      real m;
      real <lower=0> s;
      real <lower=0> a_0;
    }
    parameters{
      real mu;
    }
    model{
      /*power prior*/
        target += normal_lpdf(mu | m, s);
        target += a_0 * normal_lpdf(y0 | mu, sigma);
    }
