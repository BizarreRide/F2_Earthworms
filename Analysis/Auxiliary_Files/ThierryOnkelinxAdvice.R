library(glmmADMB)
Owls$Interaction <- interaction(Owls$FoodTreatment, Owls$SexParent)
om <- glmmadmb(SiblingNegotiation~ Interaction+(1|Nest)+offset(log(BroodSize)),zeroInflation=TRUE,family="nbinom",data=Owls)
library(multcomp)
pairwise <- glht(om, mcp(Interaction = "Tukey"))
pairwise.ci <- confint(pairwise)
library(ggplot2)
ggplot(pairwise.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + geom_errorbarh() + geom_point() + geom_vline(xintercept = 1)
# exp(estimate) = 1 :::> SiblingNegotiation = 1
par(mar=c(2,15,2,2))

plot(pairwise.ci)

head(Owls)


library(glmmADMB)
om <- glmmadmb(SiblingNegotiation~ FoodTreatment * SexParent +(1|Nest)+offset(log(BroodSize)),zeroInflation=TRUE,family="nbinom",data=Owls)
newdata <- expand.grid(
  FoodTreatment = unique(Owls$FoodTreatment),
  SexParent = unique(Owls$SexParent),
  BroodSize = 4
)
newdata <- cbind(newdata, predict(om, newdata = newdata, interval = "confidence"))

library(ggplot2)
ggplot(newdata, aes(x = FoodTreatment, colour = SexParent, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + geom_errorbar(position = position_dodge(1)) + geom_point(position = position_dodge(1))

**********************************************************************************************************
  
  ggplot(data, aes(x=age_class, y=4*anc)) +  
  geom_point(, position_dodge=0.2) +
  facet_grid(.~samcam)
  ggtitle("Data 1st Order") +
  xlab("Age Class") + 
  ylab("Abundance") +
  ylim(-10,max(data.rf$abc.mean+data.rf$abc.se)) +
  labs(fill="Functional Group") +
  scale_fill_grey(labels=c("anecic","endogeic", "epigeic", "total")) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
        legend.position=c(0.18,0.68))




