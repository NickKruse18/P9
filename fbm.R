download.file(url = "https://cran.r-project.org/src/contrib/Archive/somebm/somebm_0.1.tar.gz", destfile = "somebm_0.1.tar.gz")

install.packages(pkgs="somebm_0.1.tar.gz", type="source", repos=NULL)

library("somebm")

FBM=fbm(hurst=0.2,n=1000)
plot(FBM)
d <- fbm(hurst=0.2, n=1000)
plot(d