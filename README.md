## Lesson 7 - Tuesday 2/25/20

### Case Study #3 - Statistical Power Assessment (Continued from Last Week)

* A decision to reject a null hypothesis (or not) is a yes/no question. So, our inference will either be correct or incorrect. There are two kinds of incorrect decisions: (1) a decision to reject a hypothesis when that hypothesis is true (alpha or Type I error); and (2) a decision to accept a hypothesis when that hypothesis is false (beta or Type II error). 

* Statistical significance tests are generally calibrated in terms of keeping the probability of a Type I error small (i.e., a = 0.05 or less). 

* We would also like the probability of making a Type II error to be small. Generally, this kind of error is phrased in terms of 1 - Type II error probability which corresponds to the power of a test (i.e., the probability of rejecting a hypothesis when that hypothesis is false. 

* When designing a study it is often useful (and required) to consider the power of the test that will be conducted. Here is a standard [power calculator](https://www.stat.ubc.ca/~rollin/stats/ssize/b2.html) found on the web.

* The test we will conduct is a difference between independent two proportions z-test. Recall that all t-tests and z-tests have the following form: delta/se(delta).

* The critical value for a two-tailed z-test is abs(1.96). So, a z-test that is greater than 1.96 or less than -1.96 is in the critical region and leads us to reject the hypothesis that delta is zero in the population from which the sample came.

* Here is an example calculation of the z-test based on a single sample drawn from a large population:

```R
# set random number seed

set.seed(8479581)

# define the population

trt <- c(rep(1,50000),rep(0,50000))
out <- c(rep(0,35000),rep(1,15000),rep(0,30000),rep(1,20000))

trt.table <- table(out,trt)
trt.table

py1t1 <- trt.table[2,2]/(trt.table[1,2]+trt.table[2,2])
py1t0 <- trt.table[2,1]/(trt.table[1,1]+trt.table[2,1])
py1t1-py1t0

# create population data frame

pop.data <- data.frame(trt,out)

# draw a single sample

ncases <- 300

i <- sample(1:nrow(pop.data),size=ncases,replace=T)

sdftrt <- pop.data$trt[i]
sdfout <- pop.data$out[i]
sdftable <- table(sdfout,sdftrt)
sdftable
  
sdf.n0 <- sdftable[1,1]+sdftable[2,1]
sdf.n1 <- sdftable[1,2]+sdftable[2,2]
sdf.p0 <- sdftable[2,1]/(sdftable[1,1]+sdftable[2,1])
sdf.p1 <- sdftable[2,2]/(sdftable[1,2]+sdftable[2,2])
sdf.pooled.p <- (sdf.p0*sdf.n0+sdf.p1*sdf.n1)/(sdf.n0+sdf.n1)
sdf.delta <- sdf.p1-sdf.p0
sdf.zval <- sdf.delta/sqrt(sdf.pooled.p*(1-sdf.pooled.p)*(1/sdf.n0+1/sdf.n1)) 
sdf.delta
sdf.zval
```

Here are the results:

```Rout
> # set random number seed
> 
> set.seed(8479581)
> 
> # define the population
> 
> trt <- c(rep(1,50000),rep(0,50000))
> out <- c(rep(0,35000),rep(1,15000),rep(0,30000),rep(1,20000))
> 
> trt.table <- table(out,trt)
> trt.table
   trt
out     0     1
  0 30000 35000
  1 20000 15000
> 
> py1t1 <- trt.table[2,2]/(trt.table[1,2]+trt.table[2,2])
> py1t0 <- trt.table[2,1]/(trt.table[1,1]+trt.table[2,1])
> py1t1-py1t0
[1] -0.1
> 
> # create population data frame
> 
> pop.data <- data.frame(trt,out)
> 
> # draw a single sample
> 
> ncases <- 300
> 
> i <- sample(1:nrow(pop.data),size=ncases,replace=T)
> 
> sdftrt <- pop.data$trt[i]
> sdfout <- pop.data$out[i]
> sdftable <- table(sdfout,sdftrt)
> sdftable
      sdftrt
sdfout  0  1
     0 88 99
     1 65 48
>   
> sdf.n0 <- sdftable[1,1]+sdftable[2,1]
> sdf.n1 <- sdftable[1,2]+sdftable[2,2]
> sdf.p0 <- sdftable[2,1]/(sdftable[1,1]+sdftable[2,1])
> sdf.p1 <- sdftable[2,2]/(sdftable[1,2]+sdftable[2,2])
> sdf.pooled.p <- (sdf.p0*sdf.n0+sdf.p1*sdf.n1)/(sdf.n0+sdf.n1)
> sdf.delta <- sdf.p1-sdf.p0
> sdf.zval <- sdf.delta/sqrt(sdf.pooled.p*(1-sdf.pooled.p)*(1/sdf.n0+1/sdf.n1)) 
> sdf.delta
[1] -0.09830599
> sdf.zval
[1] -1.756649
> 
```

Note that since this z-value is greater than -1.96 and less than +1.96, it does not fall in the critical region. So, we fail to reject the null hypothesis. Since this is an example for a single sample, it is not sufficient for assessing power. We now move to a repeated sampling simulation to measure the probability of rejecting the null hypothesis when -- as in this case -- the null hypothesis is false.

```R
# begin power simulation

nsamples <- 1000

delta <- vector()
zval <- vector()

for(k in 1:nsamples){
  i <- sample(1:nrow(pop.data),size=ncases,replace=T)
  strt <- pop.data$trt[i]
  sout <- pop.data$out[i]
  stable <- table(sout,strt)
  n0 <- stable[1,1]+stable[2,1]
  n1 <- stable[1,2]+stable[2,2]
  p0 <- stable[2,1]/(stable[1,1]+stable[2,1])
  p1 <- stable[2,2]/(stable[1,2]+stable[2,2])
  pooled.p <- (p0*n0+p1*n1)/(n0+n1)
  delta[k] <- p1-p0
  zval[k] <- delta[k]/sqrt(pooled.p*(1-pooled.p)*(1/n0+1/n1))  
}

reject.zval <- ifelse(abs(zval)>1.96,1,0)
table(reject.zval)
mean(delta)
```

and here is the output:

```Rout
> # begin power simulation
> 
> nsamples <- 1000
> 
> delta <- vector()
> zval <- vector()
> 
> for(k in 1:nsamples){
+   i <- sample(1:nrow(pop.data),size=ncases,replace=T)
+   strt <- pop.data$trt[i]
+   sout <- pop.data$out[i]
+   stable <- table(sout,strt)
+   n0 <- stable[1,1]+stable[2,1]
+   n1 <- stable[1,2]+stable[2,2]
+   p0 <- stable[2,1]/(stable[1,1]+stable[2,1])
+   p1 <- stable[2,2]/(stable[1,2]+stable[2,2])
+   pooled.p <- (p0*n0+p1*n1)/(n0+n1)
+   delta[k] <- p1-p0
+   zval[k] <- delta[k]/sqrt(pooled.p*(1-pooled.p)*(1/n0+1/n1))  
+ }
> 
> reject.zval <- ifelse(abs(zval)>1.96,1,0)
> table(reject.zval)
reject.zval
  0   1 
555 445 
> mean(delta)
[1] -0.09888599
> 
```

To interpret this simulation, we create a variable called reject that is coded 1 if the null hypothesis is rejected and 0 if we fail to reject. Considering the 1000 samples we drew, we rejected the null hypothesis that delta = 0 in about 44.5% of the samples. In general, we want the rejection rate to be at least 80%. To attain that level of power, we would need to increase our sample size. Let's see what happens when we increase the number of cases from 300 to 500.

```R
# begin power simulation

nsamples <- 1000
ncases <- 500

delta <- vector()
zval <- vector()

for(k in 1:nsamples){
  i <- sample(1:nrow(pop.data),size=ncases,replace=T)
  strt <- pop.data$trt[i]
  sout <- pop.data$out[i]
  stable <- table(sout,strt)
  n0 <- stable[1,1]+stable[2,1]
  n1 <- stable[1,2]+stable[2,2]
  p0 <- stable[2,1]/(stable[1,1]+stable[2,1])
  p1 <- stable[2,2]/(stable[1,2]+stable[2,2])
  pooled.p <- (p0*n0+p1*n1)/(n0+n1)
  delta[k] <- p1-p0
  zval[k] <- delta[k]/sqrt(pooled.p*(1-pooled.p)*(1/n0+1/n1))  
}

reject.zval <- ifelse(abs(zval)>1.96,1,0)
table(reject.zval)
mean(delta)
```

When we conduct the simulation with a sample size of 500 cases instead of 300 cases, the situation improves:

```Rout
> # begin power simulation
> 
> nsamples <- 1000
> ncases <- 500
> 
> delta <- vector()
> zval <- vector()
> 
> for(k in 1:nsamples){
+   i <- sample(1:nrow(pop.data),size=ncases,replace=T)
+   strt <- pop.data$trt[i]
+   sout <- pop.data$out[i]
+   stable <- table(sout,strt)
+   n0 <- stable[1,1]+stable[2,1]
+   n1 <- stable[1,2]+stable[2,2]
+   p0 <- stable[2,1]/(stable[1,1]+stable[2,1])
+   p1 <- stable[2,2]/(stable[1,2]+stable[2,2])
+   pooled.p <- (p0*n0+p1*n1)/(n0+n1)
+   delta[k] <- p1-p0
+   zval[k] <- delta[k]/sqrt(pooled.p*(1-pooled.p)*(1/n0+1/n1))  
+ }
> 
> reject.zval <- ifelse(abs(zval)>1.96,1,0)
> table(reject.zval)
reject.zval
  0   1 
356 644 
> mean(delta)
[1] -0.09818709
> 
```

Note that we still haven't reached the threshold of acceptable power yet, but 500 cases is much closer to that level than 300 cases. For your assignment, you will need to find the minimum number of cases that yields 80% power. Confirm that your results are correct using the calculator linked above.

### Identification Issues in Criminology

#### 1. Diane Saphire (1984)

* National Crime Survey (NCS)
* p(HH not victimized in a year's time) = θnv (1975)
* S = R + NR
* p(HH participates in NCS) = R/S = p
* θnv = p x θnv|R + (1-p) x θnv|NR (Law of total probability)
* θnv|NR ∊ [0,1]
* Point-identified estimate assuming MAR = 0.732
* Interval estimate = [0.523,0.770]

#### 2. Charles Manski (1995,2012)

* Conclusions = data + assumptions
* More credible --> answer is vague
* Less credible --> answer is specific
* Law of decreasing credibility

#### Example Research Question: Did crime go up or down?

* A fairly common question: why did residential burglaries go up (or down)?
* Another, less common, question: did residential burglaries go up (or down)?
* Example: household burglaries in Charlotte NC; comparing 2014 to 2015.
* Household burglary count from CMPD in 2014 was 4,490; in 2015 the number was 4,983 (an 11% y-o-y increase).
* Some complications: population, hierarchy rule, # of people per HH, HH vacancy rate, p(report to police).
* In 2014, p(report) was estimated to be 60%; in 2015, it dropped to 50.8%.
* NCS/NCVS reports from 1973-2018 for household burglary report-to-police rates (range of 46% to 60%).
* City-specific estimates from 1970's (range of 45% to 58%). 
* Urban estimates from 1995 to 2018 (range of 46% to 62%).
* South estimates from 1996 to 2018 (range of 47% to 67%); larger y-o-y changes in last 10 years.
* Largest y-o-y change was South region from 2010-11 (64.6% to 48.4%).
* Lauritsen/Schaum estimates for NYC, Chicago, and LA (50% to 61%).
* Let p(report hh burglary to police) = pr ∊ [0.45,0.67]

```rout
# 2014 hh burglaries in CMPD

> 4490/0.67
[1] 6701
>
> 4490/0.45
[1] 9978
>

# 2015 hh burglaries in CMPD

> 4983/0.67
[1] 7437
>
> 4983/0.45
[1] 11073
> 
```

* A more common year-over-year transition would be a movement in the [0.5,0.6] range.

```rout
# 2014 hh burglaries in CMPD

> 4490/0.6
[1] 7483
>
> 4490/0.5
[1] 8980

# 2015 hh burglaries in CMPD

> 4983/0.6
[1] 8305
> 
> 4983/0.5
[1] 9966

```

* In neither instance is the sign of the y-o-y change identified; but consider this example where the p(report) range is [0.5,0.55]:

```rout
# 2014 hh burglaries in CMPD

> 4490/0.55
[1] 8164
>
> 4490/0.5
[1] 8980
>

# 2015 hh burglaries in CMPD

> 4983/0.55
[1] 9060
> 
> 4983/0.5
[1] 9966
```

Part 2 of this week's assignment: consider the residential burglaries in Charlotte-Mecklenburg from 2016 (N = 4,767) and 2017 (N = 4,240). Conduct an identification analysis for the change between these two years and assess whether we can reasonably draw any conclusions about the sign of the change from one year to the next. You can use the 2017 National Crime Victimization Survey [report](https://www.bjs.gov/content/pub/pdf/cv17.pdf) (Table 6) to assess the proportion of residential burglaries reported to the police on a national scale between the two years.
