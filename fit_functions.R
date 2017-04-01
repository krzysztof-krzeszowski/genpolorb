get_E <- function(M, e) {
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

v_get_E <- Vectorize(get_E)

fit_pd_e <- function(t, t0, p, i, e, o, p_max) {
    # Harries & Howarth 2000
    o <- o * pi / 180.0
    i <- i * pi / 180.0

    M <- (2.0 * pi / p) * (t - t0)
    M <- M %% (2 * pi)

    E <- v_get_E(M, e)
    E <- E %% ( 2 * pi)

    cos_nu <- (cos(E) - e) / (1 - e * cos(E))
    sin_nu <- (sin(E) * sqrt(1 - e ** 2)) / (1 - e * cos(E))

    sin_o_nu <- sin(o) * cos_nu + sin_nu * cos(o)

    return(p_max * (1.0 - (sin(i) ** 2) * (sin_o_nu ** 2)))
}
