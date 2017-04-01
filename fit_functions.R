get.E <- function(M, e) {
    # Most appropiate initial guesses Esmaelzadeh & Ghadiri 2014
    if (M >= 0 & M <= 0.25) {
        E0 <- M + ((6 * M) ** (1. / 3.) - M) * e
    } else if (M > 0.25 & M <= 2) {
        E0 <- M + e * (sin(M)) / (1 - sin(M + e) + sin(M))
    } else {
        E0 <- M + (e * sin(M)) / sqrt(1 - 2 * e * cos(M) + e ** 2)
    }

    E1 <- E0 - (E0 - M - e * sin(E0)) / (1 - e * cos(E0))

    i <- 0
    while (abs(E1 - E0) > 1e-10) {
        E0 <- E1
        E1 <- E0 - (E0 - M - e * sin(E0)) / (1 - e * cos(E0))
        i <- i + 1
        if (i > 10) {
            return(E1)
        }
    }
    
    return(E1)
}

v.get.E <- Vectorize(get.E)

fit.pd.e <- function(t, t0, p, i, e, o, p.max) {
    # Harries & Howarth 2000
    o <- o * pi / 180.0
    i <- i * pi / 180.0

    M <- (2.0 * pi / p) * (t - t0)
    M <- M %% (2 * pi)

    E <- v.get.E(M, e)
    E <- E %% ( 2 * pi)

    cos.nu <- (cos(E) - e) / (1 - e * cos(E))
    sin.nu <- (sin(E) * sqrt(1 - e ** 2)) / (1 - e * cos(E))

    sin.o.nu <- sin(o) * cos.nu + sin.nu * cos(o)

    return(p.max * (1.0 - (sin(i) ** 2) * (sin.o.nu ** 2)))
}
