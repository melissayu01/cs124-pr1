# CS 124 #
# PA 2 Code #

f = file.path("/Users/phillip.yu.1/Dropbox/Phillip\ Stuff/Harvard\ -\ Year\ 2/CS\ 124/Programming\ Assignments/PR2/report/results.csv")
df = read.csv(f)
attach(df)

# optimal n_0 against n
threshplot = plot(Size, Threshold, ylim = c(0, 1024),
                  ylab="optimal threshold size", 
                  xlab="size of matrix",
                main="Optimal threshold against matrix size")
abline(h=mean(df$Threshold))

# runtime against n
timeplot = plot(df$Size, df$Time, #ylim = c(0, 1024),
                ylab="average runtime (sec)", 
                xlab="size of matrix",
                main="Runtime against matrix size")
timefit = nls(Time ~ a*Size^(log(7, 2)),start=list(a=1),data=df)
summary(timefit)
new = data.frame(Size = seq(min(Size),max(Size)))
lines(new$Size,predict(timefit,newdata=new))


detach(df)
