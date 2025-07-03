### This code was used on R 4.5.0

###################################################################################
### Likelihood ratio test function for the independence hypothesis test for a
### 2x2 contingency table.
lrt.independence.test <- function(x) {
  # Computing table totals.
  n1. <- x[1,1]+x[1,2]
  n2. <- x[2,1]+x[2,2]
  n.1 <- x[1,1]+x[2,1]
  n.2 <- x[1,2]+x[2,2]
  n.. <- n1.+n2.

  # Computing the statistics of the test.
  # obs: this is already in log scale, and to avoid the "log(0)" problem,
  # the ifelse function is used. Example: 0^0 = 1, and log(1) = 0, but when
  # computing in logscale 0*log(0) = 0*"-Inf" which should be 0, but it is "NaN" in R.
  stat <- ( ifelse(   n1. == 0,0,   n1.*log(n1.))
           +ifelse(   n2. == 0,0,   n2.*log(n2.))
           +ifelse(   n.1 == 0,0,   n.1*log(n.1))
           +ifelse(   n.2 == 0,0,   n.2*log(n.2))
           -ifelse(x[1,1] == 0,0,x[1,1]*log(x[1,1]))
           -ifelse(x[1,2] == 0,0,x[1,2]*log(x[1,2]))
           -ifelse(x[2,1] == 0,0,x[2,1]*log(x[2,1]))
           -ifelse(x[2,2] == 0,0,x[2,2]*log(x[2,2]))
           -ifelse(   n.. == 0,0,   n..*log(n..)))

  # Computing the p-value using the asymptotic distribution of "stat".
  p <- pchisq(-2*stat,df=1,lower.tail=FALSE)

  return(p)
}

###################################################################################
### Likelihood ratio test function to compute the power of the
### independence hypothesis test for a 2x2 contingency table.
lrt.independence.power <- function(x,alpha=0.05,n.rep=100000) {
  # Computing table totals.
  n1. <- x[1,1]+x[1,2]
  n2. <- x[2,1]+x[2,2]
  n.1 <- x[1,1]+x[2,1]
  n.2 <- x[1,2]+x[2,2]
  n.. <- n1.+n2.
  
  # To compute the observed power, a Monte Carlo procedure is used.
  # Step 1: Compute the MLE (maximum likelihood estimates) for the model's parameters.
  theta11 <- x[1,1]/n..
  theta12 <- x[1,2]/n..
  theta21 <- x[2,1]/n..
  theta22 <- x[2,2]/n..

  count <- 0
  for (i in 1:n.rep) {
    # Step 2: Generate a 2x2 table using the MLE.
    aux <- rmultinom(1,n..,c(theta11,theta12,theta21,theta22))

    # Step 3: Count the number of times that the hypothesis was rejected.
    count <- count + ifelse(lrt.independence.test(matrix(aux,2,2,byrow=TRUE)) <= alpha,1,0)
  }

  # Step 4: The observed power of the test is approximated by the proportion of
  #         "2x2 tables" generated given the data ("under alternative hypothesis"),
  #         which rejects the independence hypothesis.
  return(count/n.rep)
}

###################################################################################
### Log-likelihood function for a Beta model.
lik.beta <- function(par,x) {
  if (all(par > 0)) {
    out <- sum(dbeta(x,shape1=par[1],shape2=par[2],log=TRUE))
  } else {
    out <- log(0)
  }
  return(out)
}

###################################################################################
### Log-likelihood function for a Beta model given the null hypothesis (H).
lik.beta.H <- function(par,x) {
  if (par > 0) {
    out <- sum(dbeta(x,shape1=par,shape2=par,log=TRUE))
  } else {
    out <- log(0)
  }
  return(out)
}

###################################################################################
### Likelihood ratio test function for the hypothesis (H) of average accuracy
### equal to 0.5, using a Beta model.
lrt.beta <- function(x) {
  # Computing the unrestricted MLE (full parametric space).
  mle <- optim(c(15,15),lik.beta,x=x,control=list(fnscale=-1),method="BFGS")
  # Computing the MLE given H (parametric space restricted by the hypothesis H).
  mle.H <- optim((mle$par[1]+mle$par[2])/2,lik.beta.H,x=x,control=list(fnscale=-1),method="Brent",lower=0,upper=1000)

  # Computing the p-value using the asymptotic distribution of the statistics of the test.
  p <- pchisq(-2*(mle.H$value-mle$value),df=1,lower.tail=FALSE)

  return(p)
}

###################################################################################
### Likelihood ratio test function to compute the power of the test for the 
### hypothesis (H) of average accuracy equal to 0.5, using a Beta model.
lrt.beta.power <- function(x,alpha=0.05,n.rep=100000) {
  # To compute the observed power, a Monte Carlo procedure is used.
  # Step 1: Computing the unrestricted MLE (full parametric space).
  mle <- optim(c(15,15),lik.beta,x=x,control=list(fnscale=-1),method="BFGS")
    
  count <- 0
  for (i in 1:n.rep) {
    # Step 2: Generate a vector from the Beta distribution using the MLE.
    aux <- rbeta(length(x),mle$par[1],mle$par[2])

    # Step 3: Count the number of times that the hypothesis was rejected.
    count <- count + ifelse(lrt.beta(aux) <= alpha,1,0)
  }

  # Step 4: The observed power of the test is approximated by the proportion of
  #         generated vectors given the data ("under alternative hypothesis"),
  #         which rejects the independence hypothesis.
  return(count/n.rep)
}


################################################
### Running

### Read data:
df3 <- read.csv("Turing3/v3_confusion_matrices.csv")
df4 <- read.csv("Turing4/v4_confusion_matrices.csv")

# Organising the data in arrays
df.array3 <- array(NA,dim=c(2,2,7))
df.array4 <- array(NA,dim=c(2,2,7))

for (i in 0:6) {
  df.array3[,,(i+1)] <- as.matrix(df3[(1+2*i):(2+2*i),2:3])
  df.array4[,,(i+1)] <- as.matrix(df4[(1+2*i):(2+2*i),2:3])
}

rownames(df.array3) <- c("Pred.H","Pred.S")
colnames(df.array3) <- c("Actual.H","Actual.S")
rownames(df.array4) <- c("Pred.H","Pred.S")
colnames(df.array4) <- c("Actual.H","Actual.S")

# set a seed for reproducibility
set.seed(1234) 

### Computing the test for each expert
res <- matrix(NA,7,6)
for (i in 1:7) {
  res[i,] <- c((df.array3[1,1,i]+df.array3[2,2,i])/sum(df.array3[,,i]),
                lrt.independence.test(df.array3[,,i]),
               lrt.independence.power(df.array3[,,i]),
               (df.array4[1,1,i]+df.array4[2,2,i])/sum(df.array4[,,i]),
                lrt.independence.test(df.array4[,,i]),
               lrt.independence.power(df.array4[,,i]))
}
print(round(res,4))

### Computing the test for the accuracy
print(round(c(lrt.beta(res[,1]),lrt.beta.power(res[,1])),3))
print(round(c(lrt.beta(res[,4]),lrt.beta.power(res[,4])),3))

