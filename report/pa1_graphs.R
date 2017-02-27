# CS 124 #
# PA 1 code #
f = file.path("/Users/phillip.yu.1/Dropbox/Phillip Stuff/Harvard - Year 2/CS 124/cs124-pr1/RandomMST/results.csv")
df = read.csv(f)
dim0=df$Dim.0

# n = powers of 2 from 7 to 17
n = c()
for (i in 7:17)
{
  n = c(n, 2^i)
}

# 0D case
dim0=df$Dim.0
dim0plot = plot(n, dim0,ylim=c(0,1.25),ylab="Average MST Weight",
                main="Dimension 0: f(n) = 1.890")
dim0fit = nls(dim0 ~ a,start=list(a=1))

summary(dim0fit)
new = data.frame(n = seq(min(n),max(n),len=200))
abline(h=1.188909)


# 2D case
dim2=df$Dim.2
dim2plot = plot(n, dim2,ylab="Average MST Weight",
                main="Dimension 2: f(n) = 0.649*sqrt(n)")
dim2fit = nls(dim2 ~ a*sqrt(n),start=list(a=1))

summary(dim2fit)
new = data.frame(n = seq(min(n),max(n),len=200))
lines(new$n,predict(dim2fit,newdata=new))

# 3D case
dim3=df$Dim.3
dim3plot = plot(n, dim3, ylab="Average MST Weight",
                main="Dimension 3: f(n) = 0.651*n^(2/3)")
dim3fit = nls(dim3 ~ a*n^(2/3),start=list(a=1))

summary(dim3fit)
new = data.frame(n = seq(min(n),max(n),len=200))
lines(new$n,predict(dim3fit,newdata=new))

# 4D case
dim4=df$Dim.4
dim4plot = plot(n, dim4, ylab="Average MST Weight",
                main="Dimension 4: f(n) = 0.689*n^(3/4)")
dim4fit = nls(dim4 ~ a*n^(3/4),start=list(a=1))

summary(dim4fit)
new = data.frame(n = seq(min(n),max(n),len=200))
lines(new$n,predict(dim4fit,newdata=new))



