annot_right <-
function(p_high, p, l, d_high, d_mid2, d_mid1, d_low, m){
    b <- max(as.matrix(unlist(p$s)))
    a <- floor((b - 1) / 4)
    spheres3d(max(m + max(p_high$f), max(p$f)) +
    p_high$s[[length(p_high$s)]] + b + 10, d_high, l,
    radius = 1, color = "black")
    text3d(max(m + max(p_high$f), max(p$f)) +
    p_high$s[[length(p_high$s)]] + 2 * b + 50, d_high, l, "= 1")
    if(d_high - d_mid2 > (1 + 1 + a)){
        spheres3d(max(m + max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + b + 10, d_mid2, l,
        radius = 1 + a, color = "black")
        text3d(max(m + max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + 2 * b + 50, d_mid2,
        l, paste(" = ", 1 + a))}
    if(d_mid2 - d_mid1 > (1 + 2 * a + 1 + a)){
        spheres3d(max(m + max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + b + 10,
        d_mid1, l, radius = 1 + 2 * a, color = "black")
        text3d(max(m + max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + 2 * b + 50, d_mid1,
        l, paste(" = ",1 + 2 * a))}
    if(d_mid1 - d_low > ((1 + 3 * a) + (1 + 2 * a))){
        spheres3d(max(m + max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + b + 10,
        d_low,l,radius = 1 + 3 * a, color = "black")
        text3d(max(m + max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + 2 * b + 50,
        d_low,l,paste(" = ", 1 + 3 * a))}
    if(d_low > ((b) + (1 + 3 * a))){
        spheres3d(max(m+max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]]+ b + 10, 0, l,
        radius = b, color = "black")
        text3d(max(m+max(p_high$f), max(p$f)) +
        p_high$s[[length(p_high$s)]] + 
        2 * b + 50, 0, l, paste(" = ", b))}
    text3d(-100, d_high, l, "Upper quartile");
    text3d(-100, d_mid2, l, "3rd quartile");
    text3d(-100, d_mid1, l, "2nd quartile");
    text3d(-100, d_low, l, "Lower quartile");
    text3d(-100, 0, l, "Overall");
    rm(list=c("a","b"))
}
