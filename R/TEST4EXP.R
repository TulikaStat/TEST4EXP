#' Gives test statistic value and critical value for one parameter exponential distribution with group size 4
#' @export
#' @param x1 is sample from Exp Population 1.
#' @param x2 is sample from Exp Population 2.
#' @param x3 is sample from Exp Population 3.
#' @param x4 is sample from Exp Population 4.
#' @param alpha is numeric value (level of significance)
TEST4EXP<-function(x1,x2,x3,x4,alpha){



        func1<-function(a,b,c) {      #ALRT
                a*log(b/c)
        }
        func2<-function(d,e){       #Heuristic
                d/e
        }

        func3<-function(f,g,h){       #LRT
                (f/h)^g
        }

        k<-4
        nsim<-rep(1,k)
        nsim[1]<-length(x1)
        nsim[2]<-length(x2)
        nsim[3]<-length(x3)
        nsim[4]<-length(x4)
        N<-sum(nsim)





        Value_ALRT1<-function(x1,x2,x3,x4){


                minimum<-0

                a1<-nsim[1]
                b1<- N*(mean(x1)-minimum)
                c1<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])+((mean(x4)-minimum)*nsim[4])

                data.frame(a1=a1,b1=b1,c1=c1)
        }
        Value_ALRT2<-function(x1,x2,x3,x4){
                minimum<-0
                a2<-nsim[2]
                b2<- N*(mean(x2)-minimum)
                c2<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])+((mean(x4)-minimum)*nsim[4])
                data.frame(a2=a2,b2=b2,c2=c2)
        }
        Value_ALRT3<-function(x1,x2,x3,x4){
                minimum<-0
                a3<-nsim[3]
                b3<- N*(mean(x3)-minimum)
                c3<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])+((mean(x4)-minimum)*nsim[4])
                data.frame(a3=a3,b3=b3,c3=c3)
        }
        Value_ALRT4<-function(x1,x2,x3,x4){
                minimum<-0
                a4<-nsim[4]
                b4<- N*(mean(x4)-minimum)
                c4<-((mean(x1)-minimum)*nsim[1])+((mean(x2)-minimum)*nsim[2])+((mean(x3)-minimum)*nsim[3])+((mean(x4)-minimum)*nsim[4])
                data.frame(a4=a4,b4=b4,c4=c4)
        }
        Value_LRT<-function(x1,x2,x3,x4){
                minimum<-0
                T_tulika<-rep(1,k)
                T_tulika[1]<-mean(x1)-minimum
                T_tulika[2]<-mean(x2)-minimum
                T_tulika[3]<-mean(x3)-minimum
                T_tulika[4]<-mean(x4)-minimum

                data.frame(lrt1= T_tulika[1],lrt2= T_tulika[2],lrt3= T_tulika[3],lrt4=T_tulika[4])
        }

        Value_Heu<-function(x1,x2,x3,x4){
                minimum<-0
                T_tulika<-rep(1,k)
                T_tulika[1]<-mean(x1)-minimum
                T_tulika[2]<-mean(x2)-minimum
                T_tulika[3]<-mean(x3)-minimum
                T_tulika[4]<-mean(x4)-minimum
                T_max<-max(T_tulika)
                T_min<-min(T_tulika)
                data.frame(T_max=T_max,T_min=T_min)
        }
        Heu_Tulika<- Value_Heu(x1,x2,x3,x4)$T_max/Value_Heu(x1,x2,x3,x4)$T_min



        Lambda_ALRT_Tulika<- -2*(func1(Value_ALRT1(x1,x2,x3,x4)$a1,Value_ALRT1(x1,x2,x3,x4)$b1,Value_ALRT1(x1,x2,x3,x4)$c1)+func1(Value_ALRT2(x1,x2,x3,x4)$a2,Value_ALRT2(x1,x2,x3,x4)$b2,Value_ALRT2(x1,x2,x3,x4)$c2)+func1(Value_ALRT3(x1,x2,x3,x4)$a3,Value_ALRT3(x1,x2,x3,x4)$b3,Value_ALRT3(x1,x2,x3,x4)$c3)+func1(Value_ALRT4(x1,x2,x3,x4)$a4,Value_ALRT4(x1,x2,x3,x4)$b4,Value_ALRT4(x1,x2,x3,x4)$c4))
        Lambda_ALRT_Tulika

        Lambda_LRT_Tulika<- N^N*( func3(Value_LRT(x1,x2,x3,x4)$lrt1, nsim[1],Value_ALRT1(x1,x2,x3,x4)$c1 )* func3(Value_LRT(x1,x2,x3,x4)$lrt2,nsim[2],Value_ALRT2(x1,x2,x3,x4)$c2 )*func3(Value_LRT(x1,x2,x3,x4)$lrt3,nsim[3],Value_ALRT3(x1,x2,x3,x4)$c3)*func3(Value_LRT(x1,x2,x3,x4)$lrt4, nsim[4],Value_ALRT4(x1,x2,x3,x4)$c4 ))
        Lambda_LRT_Tulika




        kruskal_wallis<-function(a,b){
                (12/(N*(N+1)))*(((a^2)/b))
        }


        KW<-function(x1,x2,x3,x4){
                x_bind<-c(x1,x2,x3,x4)
                rank_bind<-rank(x_bind)
                rank_x1<-sum(rank_bind[1:nsim[1]])
                rank_x2<-sum(rank_bind[nsim[1]+1:nsim[2]])
                rank_x3<-sum(rank_bind[nsim[1]+nsim[2]+1:nsim[3]])
                rank_x4<-sum(rank_bind[nsim[1]+nsim[2]+nsim[3]+1:nsim[4]])
                kw_value<-kruskal_wallis(rank_x1,nsim[1])+kruskal_wallis(rank_x2,nsim[2])+kruskal_wallis(rank_x3,nsim[3])+kruskal_wallis(rank_x4,nsim[4])-3*(N+1)

                data.frame(statistic=kw_value)

        }
        KW_Tulika<-KW(x1,x2,x3,x4)$statistic

        est_theta= ((nsim[1]*mean(x1))+(nsim[2]*mean(x2))+(nsim[3]*mean(x3))+(nsim[4]*mean(x4)))/N


        boot=1000

        Lambda_ALRT_Tulika_boot<-rep(0,boot)
        Lambda_LRT_Tulika_boot<-rep(0,boot)
        KW_Tulika_boot<-rep(0,boot)
        Heu_Tulika_boot<-rep(0,boot)
        for (j in 1:boot) {
                x1_boot<- rexp(nsim[1],1/est_theta)
                x2_boot<- rexp(nsim[2],1/est_theta)
                x3_boot<- rexp(nsim[3],1/est_theta)
                x4_boot<- rexp(nsim[4],1/est_theta)
                Lambda_ALRT_Tulika_boot[j]<- -2*(     func1(     Value_ALRT1(x1_boot,x2_boot,x3_boot,x4_boot)$a1,   Value_ALRT1(x1_boot,x2_boot,x3_boot,x4_boot)$b1,    Value_ALRT1(x1_boot,x2_boot,x3_boot,x4_boot)$c1   )+func1(    Value_ALRT2(    x1_boot,x2_boot,x3_boot,x4_boot)$a2,   Value_ALRT2(x1_boot,x2_boot,x3_boot,x4_boot)$b2,     Value_ALRT2(x1_boot,x2_boot,x3_boot,x4_boot)$c2)   +     func1(Value_ALRT3(x1_boot,x2_boot,x3_boot,x4_boot)$a3,Value_ALRT3(x1_boot,x2_boot,x3_boot,x4_boot)$b3,Value_ALRT3(x1_boot,x2_boot,x3_boot,x4_boot)$c3)   +     func1(Value_ALRT4(x1_boot,x2_boot,x3_boot,x4_boot)$a4,Value_ALRT4(x1_boot,x2_boot,x3_boot,x4_boot)$b4,Value_ALRT4(x1_boot,x2_boot,x3_boot,x4_boot)$c4) )
                Lambda_LRT_Tulika_boot[j]<- N^N*( func3(Value_LRT(x1_boot,x2_boot,x3_boot,x4_boot)$lrt1, nsim[1],Value_ALRT1(x1_boot,x2_boot,x3_boot,x4_boot)$c1 )* func3(Value_LRT(x1_boot,x2_boot,x3_boot,x4_boot)$lrt2,nsim[2],Value_ALRT2(x1_boot,x2_boot,x3_boot,x4_boot)$c2 )*func3(Value_LRT(x1_boot,x2_boot,x3_boot,x4_boot)$lrt3,nsim[3],Value_ALRT3(x1_boot,x2_boot,x3_boot,x4_boot)$c3)*func3(Value_LRT(x1_boot,x2_boot,x3_boot,x4_boot)$lrt4, nsim[4],Value_ALRT4(x1_boot,x2_boot,x3_boot,x4_boot)$c4 ))
                Heu_Tulika_boot[j]<- Value_Heu(x1_boot,x2_boot,x3_boot,x4_boot)$T_max/Value_Heu(x1_boot,x2_boot,x3_boot,x4_boot)$T_min
                KW_Tulika_boot[j]<- KW(x1_boot,x2_boot,x3_boot,x4_boot)$statistic

        }
        new_lambda_LRT<-Lambda_LRT_Tulika_boot[order(Lambda_LRT_Tulika_boot)]
        C_star_LRT<-new_lambda_LRT[floor(0.05*boot)]
        C_star_LRT
        #new_lambda_ALRT<-Lambda_ALRT_Tulika_boot[order(Lambda_ALRT_Tulika_boot)]
        #C_star_ALRT<-new_lambda_ALRT[floor(0.05*boot)]
        #C_star_ALRT
        new_lambda_Heu<-Heu_Tulika_boot[order(Heu_Tulika_boot)]
        C_star_Heu<-new_lambda_Heu[floor((1-alpha)*boot)]
        new_lambda_KW<-KW_Tulika_boot[order(KW_Tulika_boot)]
        C_star_KW<-new_lambda_KW[floor((1-alpha)*boot)]
        C_star_KW


        cv_ALRT<-qchisq(p=alpha, df=3, lower.tail=FALSE)


        data.frame(Tests=c('LRT','ALRT','U-test','KW' ), test_statistic=c(Lambda_LRT_Tulika,Lambda_ALRT_Tulika, Heu_Tulika, KW_Tulika ), critical_value=c(C_star_LRT, cv_ALRT, C_star_Heu, C_star_KW))

}
