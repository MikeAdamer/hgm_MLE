"hgm.ma.MLEso3" <-
hgm.ma.MLEso3 <- function(Y,X = c(0.01, 0.005, 0),ord=21,method="H-BFGS")
{
  ######## SETUP #########
  {
    # Creating the Pfaffian System
    P1 <- function(X)
    {
      Pfaff <- matrix(0,nrow=4,ncol=4)
      
      Pfaff[2,2] = (-X[1]*(2*X[1]^2-X[2]^2-X[3]^2))/((X[1]^2-X[3]^2)*(X[1]^2-X[2]^2))
      Pfaff[3,3] = -X[1]/(X[1]^2-X[2]^2)
      Pfaff[4,4] = -X[1]/(X[1]^2-X[3]^2)
      
      Pfaff[1,2] = 1
      Pfaff[2,3] = X[2]/(X[1]^2-X[2]^2)
      Pfaff[2,4] = X[3]/(X[1]^2-X[3]^2)
      Pfaff[3,4] = 1
      
      Pfaff[2,1] = Pfaff[1,2]
      Pfaff[3,2] = Pfaff[2,3]
      Pfaff[4,2] = Pfaff[2,4]
      Pfaff[4,3] = Pfaff[3,4]
      return(Pfaff)
    }
    
    P2 <- function(X)
    {
      Pfaff <- matrix(0,nrow=4,ncol=4)
      
      Pfaff[2,2] = X[2]/(X[1]^2-X[2]^2)
      Pfaff[3,3] = (-X[2]*(X[1]^2-2*X[2]^2+X[3]^2))/((X[2]^2-X[3]^2)*(X[1]^2-X[2]^2))
      Pfaff[4,4] = -X[2]/(X[2]^2-X[3]^2)
      
      Pfaff[1,3] = 1
      Pfaff[2,3] = -X[1]/(X[1]^2-X[2]^2)
      Pfaff[2,4] = 1
      Pfaff[3,4] = X[3]/(X[2]^2-X[3]^2)
      
      Pfaff[3,1] = Pfaff[1,3]
      Pfaff[3,2] = Pfaff[2,3]
      Pfaff[4,2] = Pfaff[2,4]
      Pfaff[4,3] = Pfaff[3,4]
      return(Pfaff)
    }
    
    P3 <- function(X)
    {
      Pfaff <- matrix(0,nrow=4,ncol=4)
      
      Pfaff[2,2] = X[3]/(X[1]^2-X[3]^2)
      Pfaff[3,3] = X[3]/(X[2]^2-X[3]^2)
      Pfaff[4,4] = (X[3]*(X[1]^2+X[2]^2-2*X[3]^2))/((X[2]^2-X[3]^2)*(X[1]^2-X[3]^2))
      
      Pfaff[1,4] = 1
      Pfaff[2,3] = 1
      Pfaff[2,4] = -X[1]/(X[1]^2-X[3]^2)
      Pfaff[3,4] = -X[2]/(X[2]^2-X[3]^2)
      
      Pfaff[4,1] = Pfaff[1,4]
      Pfaff[3,2] = Pfaff[2,3]
      Pfaff[4,2] = Pfaff[2,4]
      Pfaff[4,3] = Pfaff[3,4]
      return(Pfaff)
    }
    
    # The HGM using a Gauge transformation to avoid numerical overflows
    C_HGM_Gauge <- function(X,E)
    {
      # Create the matrix for the integration.
      RHS <- function(t, y, X)
      {
        A <- matrix(c(0,X[1],X[2],X[3],X[1],0,X[3],X[2],X[2],X[3],0,X[1],X[3],X[2],X[1],0),nrow=4)
        l_0 <- max(c(X[1]-X[2]-X[3],-X[1]+X[2]-X[3],-X[1]-X[2]+X[3],X[1]+X[2]+X[3]))
        B <- A-(2.0/t)*diag(c(0,1,1,1))-l_0*diag(c(1,1,1,1))
        return(list(B %*% y))
      }
      
      # Integrating the ODE.
      # Corollary 1 of Sei et al. is used.
      # This avoids integrating over the singular locus
      t0 <- 1e-4
      y0 <- StartPoint(t0*X,E)*exp(-t0*max(c(X[1]-X[2]-X[3],-X[1]+X[2]-X[3],-X[1]-X[2]+X[3],X[1]+X[2]+X[3])))
      sol <- ode(y0,c(t0,1),RHS,X)
      return(as.vector(sol[2,2:5]))
    }
    
    # Creating a table of the coefficients of the series expansion. 
    #This avoids costly recomputation
    E_table_var <- function(n)
    {
      n = n-1
      # From Sei et al.
      ZR <- 2
      fct <- array(0,ZR+n)
      dfct <- array(0,ZR+4*n)
      Etens <- array(0,rep(ZR+n,3))
      
      # Calculate factorials
      fct[ZR] <- 1
      for(i in 1:n){
        fct[ZR+i] <- fct[ZR+i-1]*i
      }
      # Calculate double factorials
      dfct[ZR-1] <- dfct[ZR] <- 1
      for(i in 1:(4*n)){
        dfct[ZR+i] <- dfct[ZR+i-2]*i
      }
      
      # Calculate expectation of monomials
      for(k in 0:n){
        for(l in 0:n){
          for(m in 0:n){
            if((k-l)%%2==0 && (l-m)%%2==0){
              e <- 0
              for(i in 0:(l%/%2)){ # l-n = 2i
                a0 <- choose(l,2*i)
                a1 <- dfct[ZR+k+m]*dfct[ZR+2*i-1]/dfct[ZR+k+m+2*i+1]
                a2 <- dfct[ZR+k+l-2*i-1]*dfct[ZR+2*i-1]/dfct[ZR+k+l]
                a3 <- dfct[ZR+m+l-2*i-1]*dfct[ZR+2*i-1]/dfct[ZR+m+l]
                e <- e + a0*a1*a2*a3
              }
              # Etens[ZR+k,ZR+l,ZR+m] <- e
              Etens[ZR+k-1,ZR+l-1,ZR+m-1] <- (e/fct[ZR+k]/fct[ZR+l]/fct[ZR+m])
            }else{
              Etens[ZR+k-1,ZR+l-1,ZR+m-1] <- 0
            }
          }
        }
      }
      return(Etens)
    }
    
    # Evaluating the starting point using the series expansion
    StartPoint <- function(X0=c(1,1,1),E)
    {
      max_E <- dim(E)[1]-2
      eval_c0 <- function(X0)
      {
        c0 <- 0
        for (k in 1:max_E)
        {
          for (l in 1:max_E)
          {
            for (m in 1:max_E)
            {
              c0 = c0 + X0[1]^(k-1)*X0[2]^(l-1)*X0[3]^(m-1)*E[k,l,m]
            }
          }
        }
        return(c0)
      }
      
      eval_d1c0 <- function(X0)
      {
        d1c0 <- 0
        for (k in 1:(max_E-1))
        {
          for (l in 1:max_E)
          {
            for (m in 1:max_E)
            {
              d1c0 = d1c0 + k*X0[1]^(k-1)*X0[2]^(l-1)*X0[3]^(m-1)*E[k+1,l,m]
            }
          }
        }
        return(d1c0)
      }
      
      eval_d2c0 <- function(X0)
      {
        d2c0 <- 0
        for (k in 1:max_E)
        {
          for (l in 1:(max_E-1))
          {
            for (m in 1:max_E)
            {
              d2c0 = d2c0 + l*X0[1]^(k-1)*X0[2]^(l-1)*X0[3]^(m-1)*E[k,l+1,m]
            }
          }
        }
        return(d2c0)
      }
      
      eval_d3c0 <- function(X0)
      {
        d3c0 <- 0
        for (k in 1:max_E)
        {
          for (l in 1:max_E)
          {
            for (m in 1:(max_E-1))
            {
              d3c0 = d3c0 + m*X0[1]^(k-1)*X0[2]^(l-1)*X0[3]^(m-1)*E[k,l,m+1]
            }
          }
        }
        return(d3c0)
      }
      
      
      C0 <- c(eval_c0(X0),eval_d1c0(X0),eval_d2c0(X0),eval_d3c0(X0))
      return(C0)
    }
    
    # Compute the sign preserving singular values
    SingularValues <- function(Y)
    {
      SVD <- svd(Y)
      SVD$u <- SVD$u %*% diag(c(det(SVD$u),1,1))
      SVD$v <- t(SVD$v) %*% diag(c(det(SVD$v),1,1))
      SVD$d <- diag(diag(SVD$d) %*% diag(c(det(SVD$u %*% SVD$v),1,1)))
      return(SVD)
    }
    
    # The gradient function
    Gradient <- function(X,d,E)
    {
      # First, evalutate C at a new point
      C <- C_HGM_Gauge(X,E)
      return(-d + (1/C[1])*C[2:4])
    }
    
    #  Computing the Hessian
    Hessian <- function(X,C)
    {
      # The Hessian is calculated using holonomic methods
      # Further details can be found in the reference given in the documentation
      H11 = -(P1(X) %*% C)[2]/C[1] + (1/C[1]^2)*C[2]*C[2]
      H12 = -(P2(X) %*% C)[2]/C[1] + (1/C[1]^2)*C[2]*C[3]
      H13 = -(P3(X) %*% C)[2]/C[1] + (1/C[1]^2)*C[2]*C[4]
      H21 = H12
      H22 = -(P2(X) %*% C)[3]/C[1] + (1/C[1]^2)*C[3]*C[3]
      H23 = -(P3(X) %*% C)[3]/C[1] + (1/C[1]^2)*C[3]*C[4]
      H31 = H13
      H32 = H23
      H33 = -(P3(X) %*% C)[4]/C[1] + (1/C[1]^2)*C[4]*C[4]
      Hess = matrix(c(H11,H21,H31,H12,H22,H32,H13,H23,H33),nrow = 3)
      return(Hess)
    }
    
    # The log-likelihood function which is optimised
    Likelihood <- function(X,d,E)
    {
      # First, caluculate C at a new point
      C <- C_HGM_Gauge(X,E)
      return(-d[1]*X[1]-d[2]*X[2]-d[3]*X[3] + log(C[1])+max(c(X[1]-X[2]-X[3],-X[1]+X[2]-X[3],-X[1]-X[2]+X[3],X[1]+X[2]+X[3])))
    }
    
    # storing the values of the coefficients of the series expansion
    # also compute and store the singular values of the data
    E <- E_table_var(ord)
    SVD <- SingularValues(Y)
    d <- SVD$d
  }
  
  ######## RUN #########
  {
  
  MLE <- vector(mode="list",length=2)
  names(MLE) <- c("parameter", "value")
    
  # If the singular values are too close to unity or explicitly desired the asymptotic formula is used
  if ((method=="Asymptotic") || (1-max(abs(d)) < 1e-3))
  {
    # Use the asymptotic formula
    X1 <- (0.5*(1 + 2*d[1] - 3*d[1]^2 - 2*d[2] + 2*d[1]*d[2] + d[2]^2 - 2*d[3] + 2*d[1]*d[3] - 2*d[2]*d[3] + 
          d[3]^2))/(1 - 1*d[1] - d[1]^2 + d[1]^3 - d[2] + 2*d[1]*d[2] - d[1]^2*d[2] - d[2]^2 - d[1]*d[2]^2 + 
          d[2]^3 - d[3] + 2*d[1]*d[3] - d[1]^2*d[3] + 2*d[2]*d[3] + 2*d[1]*d[2]*d[3] - d[2]^2*d[3] - 
          d[3]^2 - d[1]*d[3]^2 - d[2]*d[3]^2 + d[3]^3)
    X2 <- (0.5*(1 - 2*d[1] + d[1]^2 + 2*d[2] + 2*d[1]*d[2] - 3*d[2]^2 - 2*d[3] - 
          2*d[1]*d[3] + 2*d[2]*d[3] + d[3]^2))/(1 - d[1] - d[1]^2 + d[1]^3 - 
          d[2] + 2*d[1]*d[2] - d[1]^2*d[2] - d[2]^2 - d[1]*d[2]^2 + d[2]^3 - 
          d[3] + 2*d[1]*d[3] - d[1]^2*d[3] + 2*d[2]*d[3] + 2*d[1]*d[2]*d[3] - 
          d[2]^2*d[3] - d[3]^2 - d[1]*d[3]^2 - d[2]*d[3]^2 + d[3]^3)
    X3 <- (0.5*(1 - 2*d[1] + d[1]^2 - 2*d[2] - 2*d[1]*d[2] + d[2]^2 + 2*d[3] + 
          2*d[1]*d[3] + 2*d[2]*d[3] - 3*d[3]^2))/(1 - d[1] - d[1]^2 + d[1]^3 - 
          d[2] + 2*d[1]*d[2] - d[1]^2*d[2] - d[2]^2 - d[1]*d[2]^2 + d[2]^3 - 
          d[3] + 2*d[1]*d[3] - d[1]^2*d[3] + 2*d[2]*d[3] + 2*d[1]*d[2]*d[3] - 
          d[2]^2*d[3] - d[3]^2 - d[1]*d[3]^2 - d[2]*d[3]^2 + d[3]^3)
    
    MLE$parameter <- SVD$u %*% diag(c(X1,X2,X3)) %*% SVD$v
    MLE$value <- -Likelihood(c(X1,X2,X3),d,E)
    return(MLE)
  }
  
  # BFGS is the default procedure. The R function optim is used
  if (method=="H-BFGS")
  {
    # The tolerance is low, hence, the accuracy is limited by the numerical integration
    sol<-optim(par=X,fn=Likelihood,gr=Gradient,E=E,d=d,method=("BFGS"),control=list(reltol=1e-15,maxit=1e3))
    MLE$parameter <- SVD$u %*% diag(sol$par) %*% SVD$v
    MLE$value <- -sol$value
    return(MLE)
  }
  # Often the Newton method is the most efficient algorithm.
  else if (method=="H-Newton")
  {
    C <- StartPoint(X,E)
    g <- d - (1/C[1])*C[2:4]
    
    # Use the Newton method; the stopping condition is arbitrary and can be altered for greater accuracy.
    while(max(abs(g)) > 1e-6)
    {
      X = X - solve(Hessian(X,C))%*%g
      C = C_HGM_Gauge(X,E)
      g = d - (1/C[1])*C[2:4]
    }
    MLE$parameter <- SVD$u %*% diag(as.vector(X)) %*% SVD$v
    MLE$value <- -Likelihood(as.vector(X),d,E)
    return(MLE)
  }
  
  # In case there is a typo in the method argument
  else
  {
    warning("\'method\' must be \'Asymptotic\', \'H-BFGS\' or \'H-Newton\'.")
    return(1)
  }
  }
}
