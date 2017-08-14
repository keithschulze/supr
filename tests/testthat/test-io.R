context("Input/Output")

test_that("thunderstorm data can be read into ppp", {
    frame <- c(1, 2)
    x <- c(5177.75763, 19380.53241)
    y <- c(5767.260249, 14945.74023)
    sigma <- c(1144.627942, 105.2116346)
    intensity <- c(4426126.904, 15854.51749)
    offset <- c(3396.189696, 16694.14318)
    bkgstd <- c(1000.385184, 1484.338322)
    uncertainty <- c(8.90612651, 31.12435155)
    df <- data.frame(frame, sigma, intensity, offset, bkgstd, uncertainty)
    ppp <- spatstat::ppp(x, y, xrange = c(min(x), max(x)),
        yrange = c(min(y), max(y)), marks = df, unitname = "nm")
    ppp_read <- read.moleculelist_thunderstorm("./data/thunder.csv")

    expect_that(ppp_read, is_a("ppp"))
    expect_equal(ppp, ppp_read)
    expect_equal(colnames(spatstat::marks(ppp_read)),
        c("frame", "sigma", "intensity", "offset", "bkgstd", "uncertainty"))
})

test_that("RapidStorm data can be read into ppp", {
    x <- c(28024.4, 31148.4)
    y <- c(3214.56, 13780.9)
    frame <- c(0, 0)
    amplitude <- c(27949.7, 21886.4)
    chisq <- c(1.12972e+008, 5.5698e+007)
    bkgd <- c(7217.04, 5834.04)
    df <- data.frame(frame, amplitude, chisq, bkgd)
    ppp <- spatstat::ppp(x, y, xrange = c(min(x), max(x)),
        yrange = c(min(y), max(y)), marks = df, unitname = "nm")
    ppp_read <- read.moleculelist_rapidstorm("./data/rapid.txt")

    expect_that(ppp_read, is_a("ppp"))
    expect_equal(ppp, ppp_read)
    expect_equal(colnames(spatstat::marks(ppp_read)),
        c("frame", "amplitude", "chisq", "bkgd"))
})

test_that("multiple RapidStorm datasets can be read into ppps", {
    x <- c(28024.4, 31148.4)
    y <- c(3214.56, 13780.9)
    frame <- c(0, 0)
    amplitude <- c(27949.7, 21886.4)
    chisq <- c(1.12972e+008, 5.5698e+007)
    bkgd <- c(7217.04, 5834.04)
    df <- data.frame(frame, amplitude, chisq, bkgd)
    ppp <- spatstat::ppp(x, y, xrange=c(min(x), max(x)),
        yrange=c(min(y), max(y)), marks = df, unitname="nm")

    ppp_list <- list(ppp, ppp)
    filepaths <- c("./data/rapid.txt", "./data/rapid2.txt")
    ppp_read_list <- read.moleculelists(filepaths,
        read.moleculelist_rapidstorm)

    expect_that(ppp_read_list, is_a("list"))
    expect_equal(ppp_list, ppp_read_list)
})

test_that("New RapidStorm files with different headers can be read", {
    path <- "./data/rapid_new.txt"
    hdr <- c("x", "precision_x", "y", "precision_y", "frame", "amplitude", "chisq", "bkgd")
    x <- c(7968.46, 8161.44, 8279.78, 8738.33, 7662.00, 6666.49)
    y <- c(8205.47, 7215.19, 5469.47, 8342.27, 7521.18, 4078.94)

    ppp <- read.moleculelist_rapidstorm(path, hdrs = c("x", "precision_x", "y", "precision_y", "frame", "amplitud
e", "chisq", "bkgd"))
    df <- as.data.frame(head(ppp))
    expect_equal(df$x, x)
    expect_equal(df$y, y)
})

test_that("N-STORM dataset can be read into ppp", {
    channel.name <- c("Alexa405/Alexa647", "Alexa405/Alexa647", "Cy3/Alexa647",
        "Alexa405/Alexa647", "Cy3/Alexa647", "Non Specific Activation")
    x <- c(21837.7, 25614.3, 24397.8, 21905, 21930.5, 22805.6)
    y <- c(26818.7, 26135.6, 27213.1, 25337.6, 25343.1, 14087.2)
    xc <- c(21837.7,25614.3,24397.8,21905,21930.5, 22805.6)
    yc <- c(26818.7,26135.6,27213.1,25337.6,25343.1, 14087.2)
    height <- c(141.94905,191.6516,201.66835,123.99583,197.69766, 451.81665)
    area <- c(2640.61816,2745.26904,3570.49512,3494.45166,5489.71729, 45635.58984)
    width <- c(350.71115,352.1701,294.97519,338.3331,332.70761, 312.91064)
    phi <- c(27.76296,-23.27739,65.46878,-10.48666,88.16534, 45.87223)
    ax <- c(1.15589,1.11151,1.10812,1.06366,1.14407, 1.13600)
    bg <- c(180.14813,193.45084,192.78851,184.90379,192.74905, 200.52347)
    i <- c(2778.89063,2715.4585,3800.43701,3560.21484,5521.36914, 45796.78906)
    frame <- c(2,26,38,50,70,3)
    length <- c(3,2,5,3,5,14)
    link <- c(-1,-1,-1,-1,-1,-1)
    valid <- c(6,-3,8,7,5,-4)
    z <- c(0,0,0,0,0,0)
    zc <- c(0,0,0,0,0,0)
    photons <- c(1188.27817,1235.37107,1606.7228,1572.50325,2470.37278,20536.01543)
    lateral.localization.accuracy <- c(13.68754,11.59131,9.63432,10.03054,8.08861,1.74805)
    xw <- c(21837.7,25614.3,24397.8,21905,21930.5,22805.6)
    yw <- c(26818.7,26135.6,27213.1,25337.6,25343.1,14087.2)
    xwc <- c(21837.7,25614.3,24397.8,21905,21930.5,22805.6)
    ywc <- c(26818.7,26135.6,27213.1,25337.6,25343.1,14087.2)

    df <- data.frame(channel.name, xc, yc, height, area, width, phi, ax, bg,
        i, frame, length, link, valid, z, zc, photons,
        lateral.localization.accuracy, xw, yw, xwc, ywc)

    ppp <- spatstat::ppp(x = x, y = y, xrange = c(min(x), max(x)),
        yrange = c(min(y), max(y)), marks = df, unitname = "nm")

    ppp_read <- read.moleculelist_nstorm("./data/nstorm.txt")

    expect_that(ppp_read, is_a("ppp"))
    expect_equal(ppp_read, ppp)
    expect_equal(colnames(spatstat::marks(ppp_read)),
        c("channel.name", "xc", "yc", "height", "area", "width", "phi",
            "ax", "bg", "i", "frame", "length", "link", "valid", "z", "zc",
            "photons", "lateral.localization.accuracy", "xw", "yw",
            "xwc", "ywc"))
    expect_equal(spatstat::npoints(ppp_read), 6)

    ppp_read_ns <- read.moleculelist_nstorm("./data/nstorm.txt",
                                            nonspecific = FALSE)
    expect_equal(spatstat::npoints(ppp_read_ns), 5)
})

# test_that("Read ImageJ polygon ROIs as spatstat::win", {
# })
