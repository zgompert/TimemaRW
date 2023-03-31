## R script to test SV/fusion genotype effect on performance
load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)

## pca RW without BCTURN
pc<-prcomp(t(g_RW_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],iter.max=100,nstart=50,centers=3)
gen<-ko$cluster ## these are sorted
plot(pc$x[,1],pc$x[,2])
text(pc$x[,1],pc$x[,2],gen)


## RW source host
keep<-which(phRW$Population != 'BCTURN')
hostRW<-phRW$Host[keep]
sexRW_sub<-sexRW[which(phRW$Population != 'BCTURN')]


## lm fits, RW only 15 d weight significant
summary(lm(gemma_phRW_sub[,1] ~ gen))
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.018729   0.008554   2.190   0.0343 *
#gen         -0.008547   0.003880  -2.203   0.0333 *
#Residual standard error: 0.01932 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1058,	Adjusted R-squared:  0.08402 
#F-statistic: 4.853 on 1 and 41 DF,  p-value: 0.03328
#    Null deviance: 0.017113  on 42  degrees of freedom
#Residual deviance: 0.015302  on 41  degrees of freedom
#AIC: -213.43

summary(lm(gemma_phRW_sub[,1] ~ hostRW))
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -0.002217   0.003469  -0.639   0.5264  
#hostRWRW     0.012725   0.006858   1.855   0.0708 .
#Residual standard error: 0.01962 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.07745,	Adjusted R-squared:  0.05495 
#F-statistic: 3.442 on 1 and 41 DF,  p-value: 0.07075
#    Null deviance: 0.017113  on 42  degrees of freedom
#Residual deviance: 0.015788  on 41  degrees of freedom
#AIC: -212.09

summary(lm(gemma_phRW_sub[,1] ~ hostRW*gen))
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.006337   0.011186   0.567   0.5743
#hostRWRW      0.038318   0.020914   1.832   0.0746 .
#gen          -0.003750   0.004678  -0.802   0.4276
#hostRWRW:gen -0.019727   0.012410  -1.590   0.1200
#Residual standard error: 0.01898 on 39 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1788,	Adjusted R-squared:  0.1156
#F-statistic: 2.831 on 3 and 39 DF,  p-value: 0.0508
#    Null deviance: 0.017113  on 42  degrees of freedom
#Residual deviance: 0.014053  on 39  degrees of freedom
#AIC: -213.09

summary(lm(gemma_phRW_sub[,2] ~ gen))

#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.008068   0.011748   0.687    0.496
#gen         -0.004356   0.005329  -0.818    0.418
#Residual standard error: 0.02653 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.01604,	Adjusted R-squared:  -0.007959 
#F-statistic: 0.6683 on 1 and 41 DF,  p-value: 0.4184
#    Null deviance: 0.029331  on 42  degrees of freedom
#Residual deviance: 0.028860  on 41  degrees of freedom
#AIC: -186.15

summary(lm(gemma_phRW_sub[,2] ~ hostRW))
#             Estimate Std. Error t value Pr(>|t|)   
#(Intercept) -0.007612   0.004257  -1.788  0.08113 . 
#hostRWRW     0.026051   0.008417   3.095  0.00354 **
#Residual standard error: 0.02408 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1894,	Adjusted R-squared:  0.1696 
#F-statistic: 9.581 on 1 and 41 DF,  p-value: 0.003537
#    Null deviance: 0.029331  on 42  degrees of freedom
#Residual deviance: 0.023775  on 41  degrees of freedom
#AIC: -194.48

summary(lm(gemma_phRW_sub[,2] ~ hostRW*gen))
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.017684   0.014449  -1.224    0.228
#hostRWRW      0.038471   0.027015   1.424    0.162
#gen           0.004415   0.006042   0.731    0.469
#hostRWRW:gen -0.006029   0.016030  -0.376    0.709
#Residual standard error: 0.02452 on 39 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.2006,	Adjusted R-squared:  0.1391 
#F-statistic: 3.262 on 3 and 39 DF,  p-value: 0.03149
#    Null deviance: 0.029331  on 42  degrees of freedom
#Residual deviance: 0.023447  on 39  degrees of freedom
#AIC: -191.08


summary(lm(gemma_phRW_sub[,3] ~ gen))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.93683    0.13360   7.012 7.83e-09 ***
#gen         -0.01849    0.05999  -0.308    0.759    
#Residual standard error: 0.3088 on 47 degrees of freedom
#Multiple R-squared:  0.002017,	Adjusted R-squared:  -0.01922 
#F-statistic: 0.09499 on 1 and 47 DF,  p-value: 0.7593
#    Null deviance: 4.4898  on 48  degrees of freedom
#Residual deviance: 4.4807  on 47  degrees of freedom
#AIC: 27.846

summary(lm(gemma_phRW_sub[,3] ~ hostRW))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.89189    0.05078  17.564   <2e-16 ***
#hostRWRW     0.02477    0.10261   0.241     0.81    
#Residual standard error: 0.3089 on 47 degrees of freedom
#Multiple R-squared:  0.001239,	Adjusted R-squared:  -0.02001 
#F-statistic: 0.05829 on 1 and 47 DF,  p-value: 0.8103
#    Null deviance: 4.4898  on 48  degrees of freedom
#Residual deviance: 4.4842  on 47  degrees of freedom
#AIC: 27.885

summary(lm(gemma_phRW_sub[,3] ~ sexRW_sub))
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    0.88235    0.05285  16.697   <2e-16 ***
#sexRW_submale  0.05098    0.09551   0.534    0.596    
#Residual standard error: 0.3081 on 47 degrees of freedom
#Multiple R-squared:  0.006025,	Adjusted R-squared:  -0.01512 
#F-statistic: 0.2849 on 1 and 47 DF,  p-value: 0.596
#    Null deviance: 4.4898  on 48  degrees of freedom
#Residual deviance: 4.4627  on 47  degrees of freedom
#AIC: 27.649

summary(lm(gemma_phRW_sub[,3] ~ hostRW*gen))
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.86738    0.17816   4.869 1.42e-05 ***
#hostRWRW      0.29929    0.33645   0.890    0.378    
#gen           0.01067    0.07425   0.144    0.886    
#hostRWRW:gen -0.17734    0.19519  -0.909    0.368    
#Residual standard error: 0.3127 on 45 degrees of freedom
#Multiple R-squared:  0.02025,	Adjusted R-squared:  -0.04507 
#F-statistic:  0.31 on 3 and 45 DF,  p-value: 0.818
#    Null deviance: 4.4898  on 48  degrees of freedom
#Residual deviance: 4.3989  on 45  degrees of freedom
#AIC: 30.943

summary(lm(gemma_phRW_sub[,3] ~ sexRW_sub*gen))
#(Intercept)        0.93976    0.16001   5.873 4.82e-07 ***
#sexRW_submale     -0.02630    0.30509  -0.086    0.932    
#gen               -0.02711    0.07114  -0.381    0.705    
#sexRW_submale:gen  0.03672    0.13898   0.264    0.793    
#Residual standard error: 0.3144 on 45 degrees of freedom
#Multiple R-squared:  0.009364,	Adjusted R-squared:  -0.05668 
#F-statistic: 0.1418 on 3 and 45 DF,  p-value: 0.9344
#    Null deviance: 4.4898  on 48  degrees of freedom
#Residual deviance: 4.4478  on 45  degrees of freedom
#AIC: 31.484

summary(glm(gemma_phRW_sub[,3] ~ gen,family="binomial"))
#            Estimate Std. Error z value Pr(>|z|)  
#(Intercept)   2.6137     1.5036   1.738   0.0821 .
#gen          -0.2046     0.6520  -0.314   0.7537  

#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 32.196  on 47  degrees of freedom
#AIC: 36.196
#Number of Fisher Scoring iterations: 5

summary(glm(gemma_phRW_sub[,3] ~ hostRW,family="binomial"))
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   2.1102     0.5294   3.986 6.73e-05 ***
#hostRWRW      0.2877     1.1710   0.246    0.806    
#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 32.232  on 47  degrees of freedom
#AIC: 36.232
#Number of Fisher Scoring iterations: 5

summary(glm(gemma_phRW_sub[,3] ~ sexRW_sub,family="binomial"))
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)     2.0149     0.5323   3.785 0.000153 ***
#sexRW_submale   0.6242     1.1639   0.536 0.591791    
#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 31.978  on 47  degrees of freedom
#AIC: 35.978

summary(glm(gemma_phRW_sub[,3] ~ sexRW_sub*gen,family="binomial"))
#                  Estimate Std. Error z value Pr(>|z|)
#(Intercept)         2.5962     1.6967   1.530    0.126
#sexRW_submale      -0.2711     3.6392  -0.074    0.941
#gen                -0.2672     0.7213  -0.370    0.711
#sexRW_submale:gen   0.4214     1.6829   0.250    0.802
#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 31.828  on 45  degrees of freedom
#AIC: 39.828

summary(glm(gemma_phRW_sub[,3] ~ hostRW*gen,family="binomial"))
#              Estimate Std. Error z value Pr(>|z|)
#(Intercept)     1.8615     1.7845   1.043    0.297
#hostRWRW       33.6612  5325.7125   0.006    0.995
#gen             0.1092     0.7554   0.145    0.885
#hostRWRW:gen  -17.0659  2662.8564  -0.006    0.995
#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 30.734  on 45  degrees of freedom
#AIC: 38.734

sexRW_sub<-sexRW[which(phRW$Population != 'BCTURN')]
tapply(X=gemma_phRW_sub[,3],INDEX=list(gen,sexRW_sub),mean)
#     female  male
#1 1.0000000 1.000
#2 0.7857143 0.875
#3 0.9166667 1.000
tapply(X=gemma_phRW_sub[,3],INDEX=list(gen,sexRW_sub),sum)
#  female male
#1      8    3
#2     11    7
#3     11    4
table(list(gen,sexRW_sub))
#  female male
#  1      8    3
#  2     14    8
#  3     12    4

## pca C without BCTURN
pc<-prcomp(t(g_C_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],iter.max=100,nstart=50,centers=3)
gen<-ko$cluster ## these are sorted
plot(pc$x[,1],pc$x[,2])
text(pc$x[,1],pc$x[,2],gen)


## C source host
keep<-which(phC$Population != 'BCTURN')
hostC<-phC$Host[keep]

sexC_sub<-sexC[which(phC$Population != 'BCTURN')]
## all AIC from GLM

## lm fits, C 15 and 21 d weight and survival all significant
summary(lm(gemma_phC_sub[,1] ~ gen))
#             Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.041435   0.014971   2.768  0.00800 **
#gen         -0.018215   0.006749  -2.699  0.00958 **
#Residual standard error: 0.03677 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1317,	Adjusted R-squared:  0.1137 
#F-statistic: 7.284 on 1 and 48 DF,  p-value: 0.009578
#    Null deviance: 0.074745  on 49  degrees of freedom
#Residual deviance: 0.064898  on 48  degrees of freedom
#AIC: -184.45

summary(lm(gemma_phC_sub[,1] ~ hostC))
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0005451  0.0063400   0.086    0.932
#hostCRW     0.0125113  0.0129416   0.967    0.339
#Residual standard error: 0.03908 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.0191,	Adjusted R-squared:  -0.001336 
#F-statistic: 0.9346 on 1 and 48 DF,  p-value: 0.3385
#    Null deviance: 0.074745  on 49  degrees of freedom
#Residual deviance: 0.073318  on 48  degrees of freedom
#AIC: -178.35

summary(lm(gemma_phC_sub[,1] ~ hostC*gen))
#             Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.065608   0.023362   2.808  0.00729 **
#hostCRW     -0.055505   0.041963  -1.323  0.19247   
#gen         -0.027471   0.009538  -2.880  0.00602 **
#hostCRW:gen  0.030003   0.030018   1.000  0.32277   
#Residual standard error: 0.03674 on 46 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1691,	Adjusted R-squared:  0.1149 
#F-statistic:  3.12 on 3 and 46 DF,  p-value: 0.03493
#    Null deviance: 0.074745  on 49  degrees of freedom
#Residual deviance: 0.062107  on 46  degrees of freedom
#AIC: -182.65

summary(lm(gemma_phC_sub[,2] ~ gen))
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.045064   0.018150   2.483   0.0168 *
#gen         -0.019387   0.007982  -2.429   0.0192 *
#Residual standard error: 0.04067 on 45 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.1159,	Adjusted R-squared:  0.09625 
#F-statistic: 5.899 on 1 and 45 DF,  p-value: 0.01921
#    Null deviance: 0.084184  on 46  degrees of freedom
#Residual deviance: 0.074427  on 45  degrees of freedom
#AIC: -163.68

summary(lm(gemma_phC_sub[,2] ~ hostC))
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.002224   0.007005   0.318    0.752
#hostCRW     0.006151   0.016008   0.384    0.703
#Residual standard error: 0.04318 on 45 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.003271,	Adjusted R-squared:  -0.01888 
#F-statistic: 0.1477 on 1 and 45 DF,  p-value: 0.7026
#    Null deviance: 0.084184  on 46  degrees of freedom
#Residual deviance: 0.083908  on 45  degrees of freedom
#AIC: -158.04

summary(lm(gemma_phC_sub[,2] ~ hostC*gen))
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.05757    0.02541   2.265   0.0286 *
#hostCRW      0.03339    0.04856   0.688   0.4953  
#gen         -0.02337    0.01038  -2.252   0.0295 *
#hostCRW:gen -0.04420    0.03369  -1.312   0.1964  
#Residual standard error: 0.03997 on 43 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.1839,	Adjusted R-squared:  0.127 
#F-statistic:  3.23 on 3 and 43 DF,  p-value: 0.03153
#    Null deviance: 0.084184  on 46  degrees of freedom
#Residual deviance: 0.068701  on 43  degrees of freedom
#AIC: -163.44

summary(lm(gemma_phC_sub[,3] ~ gen))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.59451    0.12126   4.903 1.04e-05 ***
#gen          0.14099    0.05519   2.555   0.0137 *  
#---
#Residual standard error: 0.3064 on 50 degrees of freedom
#Multiple R-squared:  0.1154,	Adjusted R-squared:  0.09775 
#F-statistic: 6.526 on 1 and 50 DF,  p-value: 0.01373
#    Null deviance: 5.3077  on 51  degrees of freedom
#Residual deviance: 4.6949  on 50  degrees of freedom
#AIC: 28.522

summary(lm(gemma_phC_sub[,3] ~ hostC))
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.94872    0.04892  19.393   <2e-16 ***
#hostCRW     -0.25641    0.09784  -2.621   0.0116 *
#Residual standard error: 0.3055 on 50 degrees of freedom
#Multiple R-squared:  0.1208,	Adjusted R-squared:  0.1032
#F-statistic: 6.868 on 1 and 50 DF,  p-value: 0.01159
#    Null deviance: 5.3077  on 51  degrees of freedom
#Residual deviance: 4.6667  on 50  degrees of freedom
#AIC: 28.208

summary(lm(gemma_phC_sub[,3] ~ sexC_sub))
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.94444    0.05212   18.12   <2e-16 ***
#sexC_submale -0.19444    0.09395   -2.07   0.0437 *  
#Residual standard error: 0.3127 on 50 degrees of freedom
#Multiple R-squared:  0.0789,	Adjusted R-squared:  0.06048 
#F-statistic: 4.283 on 1 and 50 DF,  p-value: 0.04368
#    Null deviance: 5.3077  on 51  degrees of freedom
#Residual deviance: 4.8889  on 50  degrees of freedom
#AIC: 30.627

summary(lm(gemma_phC_sub[,3] ~ hostC*gen))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.83562    0.19101   4.375 6.52e-05 ***
#hostCRW     -0.56289    0.34030  -1.654    0.105    
#gen          0.04795    0.07832   0.612    0.543    
#hostCRW:gen  0.31569    0.24577   1.284    0.205    
#Residual standard error: 0.3031 on 48 degrees of freedom
#Multiple R-squared:  0.1694,	Adjusted R-squared:  0.1175 
#F-statistic: 3.264 on 3 and 48 DF,  p-value: 0.02928
#    Null deviance: 5.3077  on 51  degrees of freedom
#Residual deviance: 4.4085  on 48  degrees of freedom
#AIC: 29.248
summary(lm(gemma_phC_sub[,3] ~ sexC_sub*gen))
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       0.79240    0.15710   5.044 6.94e-06 ***
#sexC_submale     -0.43312    0.24199  -1.790   0.0798 .  
#gen               0.07018    0.06875   1.021   0.3125    
#sexC_submale:gen  0.14539    0.11545   1.259   0.2140    
#Residual standard error: 0.2997 on 48 degrees of freedom
#Multiple R-squared:  0.1879,	Adjusted R-squared:  0.1372 
#F-statistic: 3.702 on 3 and 48 DF,  p-value: 0.01779
#    Null deviance: 5.3077  on 51  degrees of freedom
#Residual deviance: 4.3103  on 48  degrees of freedom
#AIC: 28.077

summary(glm(gemma_phC_sub[,3] ~ gen,family="binomial"))

#            Estimate Std. Error z value Pr(>|z|)  
#(Intercept)  -0.8772     1.2138  -0.723   0.4699  
#gen           1.7119     0.7955   2.152   0.0314 *
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 30.603  on 50  degrees of freedom
#AIC: 34.603
#Number of Fisher Scoring iterations: 6

summary(glm(gemma_phC_sub[,3] ~ hostC,family="binomial"))
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   2.9178     0.7259   4.019 5.84e-05 ***
#hostCRW      -2.1068     0.9424  -2.236   0.0254 *  
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 31.826  on 50  degrees of freedom
#AIC: 35.826
#Number of Fisher Scoring iterations: 5

summary(glm(gemma_phC_sub[,3] ~ sexC_sub,family="binomial"))
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)    2.8332     0.7276   3.894 9.86e-05 ***
#sexC_submale  -1.7346     0.9288  -1.868   0.0618 .  
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 33.443  on 50  degrees of freedom
#AIC: 37.443


summary(glm(gemma_phC_sub[,3] ~ hostC*gen,family="binomial"))
#             Estimate Std. Error z value Pr(>|z|)
#(Intercept)    0.8769     2.4169   0.363    0.717
#hostCRW      -17.3237  2797.4433  -0.006    0.995
#gen            0.9307     1.1365   0.819    0.413
#hostCRW:gen   16.0758  2797.4422   0.006    0.995
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 29.523  on 48  degrees of freedom
#AIC: 37.523
#Number of Fisher Scoring iterations: 16
summary(glm(gemma_phC_sub[,3] ~ sexC_sub*gen,family="binomial"))
#                 Estimate Std. Error z value Pr(>|z|)
#(Intercept)        0.2023     1.9597   0.103    0.918
#sexC_submale      -1.5593     2.5392  -0.614    0.539
#gen                1.4313     1.1741   1.219    0.223
#sexC_submale:gen   0.1338     1.6052   0.083    0.934
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 28.607  on 48  degrees of freedom
#AIC: 36.607



sexC_sub<-sexC[which(phC$Population != 'BCTURN')]
tapply(X=gemma_phC_sub[,3],INDEX=list(gen,sexC_sub),mean)
#     female      male
#1 0.8571429 0.5714286
#2 0.9375000 0.8000000
#3 1.0000000 1.0000000
tapply(X=gemma_phC_sub[,3],INDEX=list(gen,sexC_sub),sum)
#female male
#1      6    4
#2     15    4
#3     13    4
table(list(gen,sexC_sub))
# female male
#  1      7    7
#  2     16    5
#  3     13    4


## combined
g_comb_sub<-cbind(g_RW_sub[a,],g_C_sub[a,])
pc<-prcomp(t(g_comb_sub),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted
save(list=ls(),file="svLm.rdat")

N_rw<-dim(gemma_phRW_sub)[1]
N_c<-dim(gemma_phC_sub)[1]
host_trt<-rep(c(1,2),c(N_rw,N_c))
summary(lm(gemma_phRW_sub[,1] ~ gen[host_trt==1]))
#                    Estimate Std. Error t value Pr(>|t|)  
#(Intercept)         0.018729   0.008554   2.190   0.0343 *
#gen[host_trt == 1] -0.008547   0.003880  -2.203   0.0333 *
#Residual standard error: 0.01932 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1058,	Adjusted R-squared:  0.08402 
#F-statistic: 4.853 on 1 and 41 DF,  p-value: 0.03328

summary(lm(gemma_phRW_sub[,2] ~ gen[host_trt==1]))
#                    Estimate Std. Error t value Pr(>|t|)
#(Intercept)         0.008068   0.011748   0.687    0.496
#gen[host_trt == 1] -0.004356   0.005329  -0.818    0.418
#Residual standard error: 0.02653 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.01604,	Adjusted R-squared:  -0.007959 
#F-statistic: 0.6683 on 1 and 41 DF,  p-value: 0.4184

summary(glm(gemma_phRW_sub[,3] ~ gen[host_trt==1],family="binomial"))
#                   Estimate Std. Error z value Pr(>|z|)  
#(Intercept)          2.6137     1.5036   1.738   0.0821 .
#gen[host_trt == 1]  -0.2046     0.6520  -0.314   0.7537  
#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 32.196  on 47  degrees of freedom
#AIC: 36.196
#Number of Fisher Scoring iterations: 5

summary(lm(gemma_phC_sub[,1] ~ gen[host_trt==2]))
#                    Estimate Std. Error t value Pr(>|t|)   
#(Intercept)         0.041435   0.014971   2.768  0.00800 **
#gen[host_trt == 2] -0.018215   0.006749  -2.699  0.00958 **
#Residual standard error: 0.03677 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1317,	Adjusted R-squared:  0.1137 
#F-statistic: 7.284 on 1 and 48 DF,  p-value: 0.009578

summary(lm(gemma_phC_sub[,2] ~ gen[host_trt==2]))
#                    Estimate Std. Error t value Pr(>|t|)  
#(Intercept)         0.045064   0.018150   2.483   0.0168 *
#gen[host_trt == 2] -0.019387   0.007982  -2.429   0.0192 *
#Residual standard error: 0.04067 on 45 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.1159,	Adjusted R-squared:  0.09625 
#F-statistic: 5.899 on 1 and 45 DF,  p-value: 0.01921

summary(glm(gemma_phC_sub[,3] ~ gen[host_trt==2],family="binomial"))
#                   Estimate Std. Error z value Pr(>|z|)  
#(Intercept)         -0.8772     1.2138  -0.723   0.4699  
#gen[host_trt == 2]   1.7119     0.7955   2.152   0.0314 *
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 30.603  on 50  degrees of freedom
#AIC: 34.603
#Number of Fisher Scoring iterations: 6

cm<-1.5;cl<-1.5

pdf("KnulliSexWeightSV.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
cs<-rep("gray10",length(sexC_sub))
a<-which(sexC_sub=="male")
cs[a]<-"red"
plot(gen[host_trt==2],gemma_phC_sub[,1],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexC_sub=="male")
o<-lm(gemma_phC_sub[a,1] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexC_sub=="female")
o<-lm(gemma_phC_sub[a,1] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="C, 15d weight",cex.main=cm)

plot(gen[host_trt==2],gemma_phC_sub[,2],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexC_sub=="male")
o<-lm(gemma_phC_sub[a,2] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexC_sub=="female")
o<-lm(gemma_phC_sub[a,2] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="C, 21d weight",cex.main=cm)


cs<-rep("gray10",length(sexRW_sub))
a<-which(sexRW_sub=="male")
cs[a]<-"red"
plot(gen[host_trt==1],gemma_phRW_sub[,1],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexRW_sub=="male")
o<-lm(gemma_phRW_sub[a,1] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexRW_sub=="female")
o<-lm(gemma_phRW_sub[a,1] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="RW, 15d weight",cex.main=cm)

plot(gen[host_trt==1],gemma_phRW_sub[,2],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexRW_sub=="male")
o<-lm(gemma_phRW_sub[a,2] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexRW_sub=="female")
o<-lm(gemma_phRW_sub[a,2] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="RW, 21d weight",cex.main=cm)
dev.off()


##################################################
csC<-c("cadetblue2","cadetblue3","cadetblue4")
csRW<-c("firebrick2","firebrick3","firebrick4")
cl<-1.6;lx<-2.1;cm<-1.6;ca<-1.2;cx<-1.3

pdf("knulliPerformanceSV.pdf",width=12,height=12)
layout(matrix(c(1,1,1,2:7),nrow=3,ncol=3,byrow=TRUE),widths=c(4,4,4),height=c(4,4,4))

plot(c(0,1),type='n',xlab="",ylab="",axes=FALSE)
title("(a) Experimental design",cex.main=cm)



par(mar=c(4.5,5.5,2.5,1.5))
plot(jitter(gen[host_trt==2]-1),gemma_phC_sub[,1],pch=19,col=csC[gen[host_trt==2]],axes=FALSE,xlab="Genotype",ylab="Resid. 15-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phC_sub[,1],INDEX=gen[host_trt==2],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.031",side=3,adj=.1,line=-2)
title("(b) 15d weight, C",cex.main=cm)
box()

plot(jitter(gen[host_trt==2]-1),gemma_phC_sub[,2],pch=19,col=csC[gen[host_trt==2]],axes=FALSE,xlab="Genotype",ylab="Resid. 21-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phC_sub[,2],INDEX=gen[host_trt==2],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.002",side=3,adj=.1,line=-2)
title("(c) 21d weight, C",cex.main=cm)
box()

surv<-tapply(INDEX=gen[host_trt==2]-1,gemma_phC_sub[,3],mean)
barplot(surv,col=csC,xlab="Genotype",ylab="Survival proportion",cex.lab=cl,cex.axis=ca,cex=cx)
mtext("P = 0.031",side=3,adj=.5,line=-2)
title("(d) Survival, C",cex.main=cm)

plot(jitter(gen[host_trt==1]-1),gemma_phRW_sub[,1],pch=19,col=csRW[gen[host_trt==1]],axes=FALSE,xlab="Genotype",ylab="Resid. 15-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phRW_sub[,1],INDEX=gen[host_trt==1],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.034",side=3,adj=.1,line=-2)
title("(e) 15d weight, RW",cex.main=cm)
box()


plot(jitter(gen[host_trt==1]-1),gemma_phRW_sub[,2],pch=19,col=csRW[gen[host_trt==1]],axes=FALSE,xlab="Genotype",ylab="Resid. 21-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phRW_sub[,2],INDEX=gen[host_trt==1],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.418",side=3,adj=.1,line=-2)
title("(f) 21d weight, RW",cex.main=cm)
box()


surv<-tapply(INDEX=gen[host_trt==1]-1,gemma_phRW_sub[,3],mean)
barplot(surv,col=csRW,xlab="Genotype",ylab="Survival proportion",cex.lab=cl,cex.axis=ca,cex=cx)
mtext("P = 0.754",side=3,adj=.5,line=-2)
title("(g) Survival, RW",cex.main=cm)
dev.off()


###################################################
dat<-dat[1:138,]
host_source_RW<-as.character(dat$Host[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "RW" & dat$Species == "knulli")])
host_source_C<-as.character(dat$Host[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "C" & dat$Species == "knulli")])
host_source<-c(host_source_RW,host_source_C)

dat_sub_C<-dat[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "C" & dat$Species == "knulli"),]
dat_sub_RW<-dat[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "RW" & dat$Species == "knulli"),]


###################################################

tapply(X=gen==1,INDEX=host_source,mean)
#        C        RW 
#0.1052632 0.6800000 
tapply(X=gen==2,INDEX=host_source,mean)
#        C        RW 
#0.4605263 0.3200000 
tapply(X=gen==3,INDEX=host_source,mean)
#        C        RW 
#0.4342105 0.0000000 

