library(MDFS)

X<- madelon$data

y<- madelon$decision

Xb<- apply(X,2, function(x) (x > median(x))*1. )

constant_cols<- apply(Xb, 2, function(x) all(x== x[[1]]) )

Xb<- Xb[, !constant_cols ]

minority_counts<-apply(Xb,2, function(x) min(table(x)))

Xb<- Xb[ , minority_counts > 30 ]

# IG2D( Z, Y)= max_j H(Y| X_j) - H(Y| X_j, Z)
# IG1D(X_j | Z)= H( Y| Z) - H( Y| X_j,Z)
# H( Y| X_j,Z) =H( Y| Z) - IG1D(X_j | Z)

#IG1D(Z)= H(Y) - H(Y|Z)
#H(Y|Z) = H(Y) - IG1D(Z)
# so also H(Y| X_j) = H(Y) - IG1D(X_j)

#therefore
# IG2D(Z,Y)= max_j{ H(Y) - IG1D(X_j)  - H(Y|Z) + IG1D(X_j | Z) } 
# IG2D(Z,Y)= max_j{IG1D(Z) - IG1D(X_j) + IG1D(X_j | Z) }
# DELTA IG(Z):= IG2D(Z) - IG1D(Z) = max_j { IG1D(X_j | Z) - IG1D(X_j) }
MDFS.discrete(data=Xb, decision=y,dimensions = 2, 
    pc.xi = 0.25,seed=123)-> mdfs_2d 

IG2D<- mdfs_2d$statistic
IG2D_rel<- IG2D[ mdfs_2d$relevant.variables ]

MDFS.discrete(data=Xb, decision=y,dimensions = 1, 
    pc.xi = 0.25,seed=123)-> mdfs_1d 

IG1D<- mdfs_1d$statistic


for (i in 1:5) {
IG2D_Z<- IG2D_rel[[i]]



IG1D_Z<- IG1D[ mdfs_2d$relevant.variables ][[i]]

DELTA_Z<- IG2D_Z - IG1D_Z

IG1D_givenZ<- MDFS.discrete.confounded(data= Xb, decision=y, confounders= Xb[, mdfs_2d$relevant.variables[[i]] ], dimensions=1)
IG1D_givenZ<- IG1D_givenZ$statistic
is_it_delta_Z<-max( IG1D_givenZ - IG1D )
print(sprintf("variable %d",i))
print(sprintf("delta_IG (IG2D - IG1D)= %.6f", DELTA_Z))
print(sprintf("delta_IG (max_j  IG1D(X_j | Z) - IG1D(X_j) ) = %.6f", is_it_delta_Z))
}




    


