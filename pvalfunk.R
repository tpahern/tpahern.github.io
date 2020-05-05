#Function for calculating and graphing confidence interval functions
 #written by Tom Ahern, May 2016 (02tahern@med.uvm.edu)
 #required package: ggplot2

cifunk <- function(rrpoint,lcl95,ucl95,minhyp,maxhyp,plotpoints){
  
  #derive standard error from confidence limits
  selnrr <- (log(ucl95)-log(lcl95))/3.92
  
  #create vector of hypothesized relative risks
  hyps_ <- seq(log(minhyp),log(maxhyp),length=plotpoints)
  hyps <- c(hyps_,log(rrpoint))
  hyprrs <- exp(hyps)
  
  #calculate z-scores for point-hyps
  zvals <- (log(rrpoint)-hyps)/selnrr
  
  #calculate 2-sided p-values
  pvals <- 2*(1-pnorm(abs(zvals),0,1))
  
  #combine hypotheses and 2-sided pvals in data frame
  plotinput <- data.frame(hyprrs,pvals)
  print("Plot input")
  print(plotinput)
  
  #make plot of confidence interval function
  require(ggplot2)
  plot <- ggplot(plotinput, aes(x=hyprrs, y=pvals)) + geom_line() + geom_point() + 
    geom_hline(aes(yintercept=0.05), color="lightsteelblue") +
    annotate("text", min(plotinput$hyprrs)+0.1, 0.1, label = "95% CI", color="lightsteelblue") +
    geom_vline(aes(xintercept=1.00), color="darkgrey") +
    scale_x_continuous(name="Hypothesized relative risk") +
    scale_y_continuous(name="2-sided p-value")
  return(plot)
}