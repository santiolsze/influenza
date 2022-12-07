library(dplyr)
library(tidyverse)
library(tidymodels)
library(tidyr)
library(modelr)
library(car)
library(cowplot)
library(pROC)
library(tableone)

rm(list = ls())
gc(verbose = FALSE)
setwd("/home/santiago/Documents/MinSalNacion/influenza/bases_segundo_informe/")
datapath = "flu_ve_2022_11_10_11_11_49_.csv"
dataname = "Controles_todo_negativo"

odds_as_df <- function(model){
  coef <- model$coefficients
  confs <- confint(model)
  data.frame( status = rownames(confs),
              Odds = round(exp(coef),3), 
              CILow = round(exp(confs[,"2.5 %"]),3), 
              CIHigh = round(exp(confs[,"97.5 %"]),3))}

get_n_cases <- function(dataset){
  dataset[,"status"] %>% sum()
}

deviance_explicada <-function(model){
  (model$null.deviance - model$deviance) / (model$null.deviance ) * 100
}
Hosmer_Lemeshow_plot <- function(dataset, predicted_column, class_column, bins, positive_value, color='forestgreen', nudge_x=0, nudge_y=0.05){
  # Asignar los grupos a las observaciones de acuerdo a la probabilidad predicha
  dataset['group'] <- bin(dataset[predicted_column], nbins = bins, method = 'l', labels=c(1:bins))
  # Contar la cantidad de casos positivos por grupo
  positive_class <- dataset %>% filter(!!sym(class_column)==positive_value) %>% group_by(group) %>% count()
  # Obtener la media de las predicciones por grupo
  HL_df <- dataset %>% group_by(group) %>% summarise(pred=mean(!!sym(predicted_column)), count=n()) %>%
    inner_join(.,positive_class) %>%
    mutate(freq=n/count)
  # Gráfico 
  HM_plot <- ggplot(HL_df, aes(x=pred, y=freq)) + 
    geom_point(aes(size=n), color=color) +
    geom_text(aes(label=n),nudge_y = nudge_y)+
    geom_abline(slope = 1, intercept = 0, linetype='dashed') + 
    theme_bw() +
    labs(title='Hosmer-Lemeshow', size='Casos', x="Probabilidad Predicha", y="Frecuencia observada")
  return(HM_plot)
}



fludata <- read.csv(datapath)

# Defino grupos etarios y aprovecho para eliminar a los >100.
fludata <-fludata %>% mutate(
  
  GRUPO_EDAD = case_when(
    EDAD_DIAGNOSTICO <= 2 ~ '0-2',
    EDAD_DIAGNOSTICO <= 18 ~ '03-18',
    EDAD_DIAGNOSTICO <= 50 ~ '19-50',
    EDAD_DIAGNOSTICO <= 75 ~ '51-75',
    EDAD_DIAGNOSTICO <= 100 ~ '76-100',
    T ~ '>100'
    
  )
  
) %>% filter(GRUPO_EDAD != '>100')



logit_formulas <- formulas(.response = ~ status,
                           fluvac = ~ fluvac,
                           fluvac_mes = ~ fluvac + MES_AÑO_ESTUDIO,
                           fluvac_grupoedad_sexo = ~fluvac + GRUPO_EDAD + SEXO,
                           fluvac_edad_sexo = ~fluvac + EDAD_DIAGNOSTICO + SEXO,
                           fluvac_grupoedad_sexo_mes = ~fluvac + GRUPO_EDAD + SEXO + MES_AÑO_ESTUDIO,
                           #fluvac_grupoedad_sexo_mes_int = ~fluvac*GRUPO_EDAD + SEXO + MES_AÑO_ESTUDIO,
                           fluvac_edad_sexo_mes = ~fluvac + EDAD_DIAGNOSTICO + SEXO + MES_AÑO_ESTUDIO,
                           )


models <- data_frame(logit_formulas) %>% # dataframe a partir del objeto formulas
  mutate(models = names(logit_formulas), # columna con los nombres de las formulas
         expression = paste(logit_formulas), # columna con las expresiones de las formulas
         mod = map(logit_formulas, ~glm(., data =fludata, family = "binomial"))) 

models <- models %>% 
  mutate(glance = map(mod,glance))

# Acá veo según AIC, BIC y DevExp que el mejor es exclusivamente uno.
models %>%  unnest(glance) %>% mutate(DEVEXP = (null.deviance - deviance)/null.deviance) %>% select(expression,AIC,BIC,DEVEXP) %>% 
  arrange(AIC)

# Puedo ver que no hay VIF elevado
vif(models$mod$fluvac_grupoedad_sexo_mes)

# Veo las estimaciones
odds_as_df(models$mod$fluvac_grupoedad_sexo_mes)

# Agrego las predicciones
models <- models %>% 
  mutate(pred= map(mod, augment, type.predict = "response", data = fludata))


final <- models %>% 
  filter(models=="fluvac_grupoedad_sexo_mes") %>% 
  unnest(pred)


roc_full <- roc(response=final$status, predictor=final$.fitted)

viol <- ggplot(final, aes(x=status, y=.fitted, group=status, fill=factor(status))) + 
  geom_violin() +
  theme_bw() +
  guides(scale="none", fill=guide_legend(title="Caso"))+
  labs(title='Violin plot', subtitle='fluvac_grupoedad_sexo_mes', y='Predicted probability')

hl <- Hosmer_Lemeshow_plot(final, '.fitted', 'status', 20, 1) +
    labs(subtitle="fluvac_grupoedad_sexo_mes") + xlim(0,1) + ylim(0,1) + theme(legend.position =  "None")

roc_plot <- ggroc(list(full=roc_full), size=1) + 
  geom_abline(slope = 1, intercept = 1, linetype='dashed') +
  theme_bw() + 
  labs(title='ROC curve',subtitle = paste("AUC:",round(roc_full$auc,3))) + theme(legend.position =  "None")

plot_grid(viol, hl,roc_plot, nrow = 3)


################### Subgroup modelling #########################################

fludatas <- list("0-100" = fludata,
                 #"0 años" = fludata %>% filter(EDAD_DIAGNOSTICO == '0'),
                #"1-2 años" = fludata %>% filter(EDAD_DIAGNOSTICO %in% c(1,2)),
                "0-2" = fludata %>% filter(GRUPO_EDAD == '0-2'),
                 "03-18" = fludata %>% filter(GRUPO_EDAD == '03-18'),
                 "19-50" = fludata %>% filter(GRUPO_EDAD == '19-50'),
                 "51-75" = fludata %>% filter(GRUPO_EDAD == '51-75'),
                 "76-100" = fludata %>% filter(GRUPO_EDAD == '76-100')
                #, "Influenza A" =  fludata %>% filter((grepl("Influenza A",CLASIFICACION_ALGORITMO)|(grepl("sin resultados positivos",CLASIFICACION_ALGORITMO)))),
                #  "Influenza B" = fludata %>% filter((grepl("Influenza B",CLASIFICACION_ALGORITMO)|(grepl("sin resultados positivos",CLASIFICACION_ALGORITMO))))
                 )

subgroup_models <- data_frame(fludatas) %>% 
  mutate(name = names(fludatas),
         n = map(fludatas,nrow(.)),
         mod = map(fludatas, ~glm(status~fluvac  + SEXO + MES_AÑO_ESTUDIO, data =., family = "binomial")))


subgroup_models <- subgroup_models %>% 
  mutate(glance = map(mod,glance),
         odds_as_df = map(mod, odds_as_df))

# Acá veo según AIC, BIC y DevExp que el mejor es exclusivamente uno.
subgroup_models %>%  unnest(glance) %>% mutate(DEVEXP = (null.deviance - deviance)/null.deviance) %>% select(name,AIC,BIC,DEVEXP) %>% 
  arrange(AIC)

# Puedo ver que no hay VIF elevado
vif(subgroup_models$mod$`19-50`)

# Veo las estimaciones

subgroup_models %>% 
  unnest(odds_as_df) %>%
  filter(status == "fluvacVac2022") %>% 
  select(name,n, Odds,CILow, CIHigh) %>% 
  write_csv(paste0(dataname,"_","results.csv"))

subgroup_models %>% 
  unnest(odds_as_df) %>%
  filter(status == "fluvacVacPre2022") %>% 
  select(name,n, Odds,CILow, CIHigh) %>% 
  write_csv(paste0(dataname,"_","results_pre2022.csv"))



t <- CreateTableOne(strata = c("status"),
                    vars = c("GRUPO_EDAD","SEXO","fluvac","MES_AÑO_ESTUDIO","CLASIFICACION_ALGORITMO"),
                     factorVars = c("GRUPO_EDAD","SEXO","fluvac","MES_AÑO_ESTUDIO","CLASIFICACION_ALGORITMO")
                  , data = fludata)

t
write_csv(as.data.frame(print(t, showAllLevels = TRUE)),paste0(dataname,"_","descriptive.csv"))


for (subdf in fludatas){
  u <- unique(subdf$GRUPO_EDAD)
  if (length(u) == 1){
    name <- paste0("raw_or_",unique(subdf$GRUPO_EDAD),dataname,".csv")  
  }else{name <- paste0("raw_or_","0-100",dataname,".csv")}
  
  tabla <- table(subdf$fluvac,subdf$status)
  as.data.frame(oddsratio(tabla)$measure) %>% mutate(level = rownames(.)) %>% write_csv(name)
}


library(ggpmisc)

fludata %>% 
  group_by(EDAD_DIAGNOSTICO) %>%
  count(status) %>%
  mutate(prop = n/sum(n)) %>%
  filter(status == 1) %>%
  summarise(log_odds = log(prop/(1 - prop))) %>%
  ggplot(aes(x = EDAD_DIAGNOSTICO, y = log_odds)) +
  geom_point() +geom_smooth(method = "lm") + scale_x_continuous(n.breaks = 30) +
  ylab("Log odds of hosp") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +theme_bw() +
  stat_poly_line() +
  stat_poly_eq()


