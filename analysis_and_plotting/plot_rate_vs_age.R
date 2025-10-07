library(ggplot2)

df = read.csv("age_effects.tsv", sep="\t")

g <- facet_wrap()


m = glm(data=subset(df, POI == "dad" & TR == "non-homopolymer STR"), count ~ PaAge, family=poisson(link="identity"))
print (summary(m))