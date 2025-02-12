# true loading matrices with 3, 5, and 7 factors
Ptrue3_f <- matrix(c(.7, .8, .9, rep(0,9), .7, .8, .9, rep(0,9), .7, .8, .9), nrow = 9, ncol = 3, byrow = F)
Ptrue5_f <- matrix(c(.6, .6, .7, .8, .8, rep(0, 15), .6, .6, .7, .8, .8, rep(0, 15), .6, .6, .7, .8, .8), ncol = 3)
Ptrue7_f <- matrix(c(.6, .6, .7, .7, .7, .8, .8, rep(0,21), .6, .6, .7, .7, .7, .8, .8, rep(0,21), .6, .6, .7, .7, .7, .8, .8), ncol = 3)

# Correlation matrix with 3 factors
corr3_f <- matrix(c(1, .7*.8, .7*.9, rep(0,6), -.02*.7,
                    .7*.8, 1, .8*.9, rep(0,6), -.02*.8,
                    .7*.9, .9*.8, 1, rep(0,6), -.02*.9,
                    rep(0,3), 1, .7*.8, .7*.9, rep(0,3), .15*.7,
                    rep(0,3), .7*.8, 1, .8*.9, rep(0,3), .15*.8,
                    rep(0,3), .7*.9, .9*.8, 1, rep(0,3), .15*.9,
                    rep(0,6), 1, .7*.8, .7*.9, .7*.5,
                    rep(0,6), .7*.8, 1, .8*.9, .8*.5,
                    rep(0,6), .7*.9, .9*.8, 1, .9*.5,
                    -.02*.7, -.02*.8, -.02*.9, .15*.7, .15*.8, .15*.9, .5*.7, .5*.8, .5*.9, 1
), 10, 10, byrow = T)

# Correlation matrix with 5 factors
corr5_f <- matrix(c(1, .6^2, .6*.7, .6*.8, .6*.8, rep(0,10), -.02*.6,
                    .6^2, 1, .6*.7, .6*.8, .6*.8, rep(0,10), -.02*.6,
                    .6*.7, .6^2, 1, .7*.8, .7*.8, rep(0,10), -.02*.7,
                    .8*.6, .8*.6, .8*.7, 1, .8*.8, rep(0,10), -.02*.8,
                    .8*.6, .8*.6, .8*.7, .8*.8, 1, rep(0,10), -.02*.8,
                    rep(0,5), 1, .6^2, .6*.7, .6*.8, .6*.8, rep(0,5), .15*.6,
                    rep(0,5), .6^2, 1, .6*.7, .6*.8, .6*.8, rep(0,5), .15*.6,
                    rep(0,5), .6*.7, .6^2, 1, .7*.8, .7*.8, rep(0,5), .15*.7,
                    rep(0,5), .8*.6, .8*.6, .8*.7, 1, .8*.8, rep(0,5), .15*.8,
                    rep(0,5), .8*.6, .8*.6, .8*.7, .8*.8, 1, rep(0,5), .15*.8,
                    rep(0,10), 1, .6^2, .6*.7, .6*.8, .6*.8, .5*.6,
                    rep(0,10), .6^2, 1, .6*.7, .6*.8, .6*.8, .5*.6,
                    rep(0,10), .6*.7, .6^2, 1, .7*.8, .7*.8, .5*.7,
                    rep(0,10), .8*.6, .8*.6, .8*.7, 1, .8*.8, .5*.8,
                    rep(0,10), .8*.6, .8*.6, .8*.7, .8*.8, 1, .5*.8,
                    -.02*.6, -.02*.6, -.02*.7, -.02*.8, -.02*.8,
                    .15*.6, .15*.6, .15*.7, .15*.8, .15*.8,
                    .5*.6, .5*.6, .5*.7, .5*.8, .5*.8, 1), 16, 16, byrow = T)

# Correlation matrix with 7 factors

PSI1 <- diag(1-rowSums(Ptrue7f^2))
SIGMA1 <- Ptrue7f%*%t(Ptrue7f)+PSI1 
y <- c(-.02*.6, -.02*.6, -.02*.7, -.02*.7, -.02*.7, -.02*.8, -.02*.8,
       .15*.6, .15*.6, .15*.7, .15*.7, .15*.7, .15*.8, .15*.8,
       .5*.6, .5*.6, .5*.7, .5*.7, .5*.7, .5*.8, .5*.8)
ty <- c(-.02*.6, -.02*.6, -.02*.7, -.02*.7, -.02*.7, -.02*.8, -.02*.8,
        .15*.6, .15*.6, .15*.7, .15*.7, .15*.7, .15*.8, .15*.8,
        .5*.6, .5*.6, .5*.7, .5*.7, .5*.7, .5*.8, .5*.8, 1)
a <- cbind(SIGMA1, c(-.02*.6, -.02*.6, -.02*.7, -.02*.7, -.02*.7, -.02*.8, -.02*.8,
                     .15*.6, .15*.6, .15*.7, .15*.7, .15*.7, .15*.8, .15*.8,
                     .5*.6, .5*.6, .5*.7, .5*.7, .5*.7, .5*.8, .5*.8))
corr7_f <- rbind(a, c(-.02*.6, -.02*.6, -.02*.7, -.02*.7, -.02*.7, -.02*.8, -.02*.8,
                      .15*.6, .15*.6, .15*.7, .15*.7, .15*.7, .15*.8, .15*.8,
                      .5*.6, .5*.6, .5*.7, .5*.7, .5*.7, .5*.8, .5*.8, 1))





