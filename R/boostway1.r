boostway1 = function(x, y, weight, nlambda, lambda.min, standardize, epsilon)
{

nobs=nrow(x)
np=ncol(x)

x.flag=0
x.name = colnames(x)
if (is.null(x.name)) {x.flag=1;x.name = paste("V", seq(np), sep = "")}

x = as.matrix(x)
if (standardize) 
	{ 
	x.train = scale(x) 
	center = attributes(x.train)$`scaled:center`
	scale = attributes(x.train)$`scaled:scale`
	} else { x.train=x }

y.flag=0
if (is.numeric(y)) {y.flag=1}
y.name = levels(factor(y))
k = length(y.name)
if ( as.integer(k)<2 ) stop("y should be of more than 1 class")
kminus = k-1
y.temp = as.factor(y)
y.train = rep(0,length(y))	
for (ii in 1:k){ y.train[which(y.temp %in% y.name[ii])]=ii }

Y.matrix = Y.matrix.gen(y.train,k)

if (missing(weight)) { weight = rep(1, nobs) }

warmbeta = matrix(rep(0,((np+1)*kminus)),(np+1),kminus)
warminner = rep(0,nobs)

aa = .C("anglebooststart",as.vector((x.train)),as.vector((Y.matrix)),
	as.integer(kminus),as.integer(nobs),as.integer(np),
	as.double(epsilon),as.double(weight),as.vector((warmbeta)),
	as.double(warminner),
	betaout=as.vector(matrix(rep(0,((np+1)*kminus)),(np+1),kminus)),
	innerout=as.double(rep(0,nobs)),PACKAGE='smac'
	)

warmbeta = aa$betaout
warminner = aa$innerout

tempbetaout=matrix(warmbeta,(np+1),(k-1))
tempbeta0 = tempbetaout[1,]
tempbeta = tempbetaout[-1,]

zerofd = matrix(0, np, kminus)
tempfd = firboost(warminner)
	for (p in 1:np)
	{
		for (q in 1:kminus)
		{
		zerofd[p,q] = sum(tempfd*x.train[,p]*Y.matrix[,q]*weight)/nobs
		}
	}

if (standardize) {
tempbeta0 = as.vector(tempbeta0 - t(as.matrix(center/scale))%*%tempbeta)
tempbeta = tempbeta/scale
}

if (inherits(tempbeta, "numeric")) {tempbeta=(as.matrix(tempbeta))}

rownames(tempbeta)=x.name

lambda.trace = max(abs(zerofd))
beta0.trace = list(tempbeta0)
beta.trace = list(tempbeta)

#---------------------------------------------#

if (as.integer(nlambda)>1)
{

lambda.panel = exp(seq(log(lambda.trace),log(lambda.trace*lambda.min),len=nlambda))[-1]

for (zz in 1:length(lambda.panel))
{
aa = .C("angleboost",as.vector((x.train)),as.vector((Y.matrix)),
	as.integer(kminus),as.integer(nobs),as.integer(np), 
	as.double(lambda.panel[zz]),
	as.double(epsilon),as.double(weight),as.vector((warmbeta)),
	as.double(warminner),
	betaout=as.vector(matrix(rep(0,((np+1)*kminus)),(np+1),kminus)),
	innerout=as.double(rep(0,nobs)),PACKAGE='smac'
	)

warmbeta = aa$betaout
warminner = aa$innerout

tempbetaout=matrix(warmbeta,(np+1),(k-1))
tempbeta0 = tempbetaout[1,]
tempbeta = tempbetaout[-1,]

if ( is.na(sum(aa$betaout)) | is.nan(sum(aa$betaout)) ) {break}

if (standardize) {
tempbeta0 = as.vector(tempbeta0 - t(as.matrix(center/scale))%*%tempbeta)
tempbeta = tempbeta/scale
}

if (inherits(tempbeta, "numeric")) {tempbeta=(as.matrix(tempbeta))}

rownames(tempbeta)=x.name

lambda.trace = c(lambda.trace, lambda.panel[zz])
beta0.trace = c(beta0.trace, list(tempbeta0))
beta.trace = c(beta.trace, list(tempbeta))
}

}

if (y.flag==1) {y.name=as.numeric(y.name)}

if (x.flag) {x=data.frame(x);colnames(x)=x.name}

z=list(x=x,y=y,weight=weight,standardize=standardize,epsilon=epsilon,
k=k,x.name = x.name, y.name = y.name, lambda=lambda.trace,
beta0=beta0.trace,beta=beta.trace,loss="boost",way=1)
#---------------------------------------------#

return(z)

}
