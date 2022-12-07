rm(list = ls())
gc(verbose = FALSE)
library(dplyr)
library(lubridate)
library(car) 
library(epitools)

dates_to_numeric <- function(dataset, date_cols){
  for (date_col in date_cols){
    dataset[,date_col] <- as.numeric(unlist(dataset[,date_col]))
  }
  return(dataset)
}
################# PARÁMETROS ##################################
#"SARS COV-2 Negativo por test de antígeno"                    "Estudiados por pruebas moleculares sin resultados positivos"
#"Virus sincicial respiratorio (VSR)"                          "Sin clasificar"                                             
#"Otro virus respiratorio"                                     "SARS-COV-2 por métodos moleculares"                         
#"Co-detección de SARS COV 2 y OVR"                            "Co-detección de VSR y OVR"                                  
#"Influenza B - linaje Victoria"                               "SARS-COV-2 por test de antígeno"                            
# "Influenza A - sin subtipificar"                              "Influenza B - sin linaje"                                   
# "Influenza A - H3N2"                                          "Co-detección de Influenza y SARS COV 2"                     
# "Co-detección de Influenza y OVR"                             "Co-detección de SARS COV 2 y VSR"                           
# "Influenza A - H1N1"                                          "Co-detección de Influenza y VSR"                            
# "Caso invalidado por epidemiología"            

clasif_casos <- c("Influenza A - sin subtipificar", "Influenza A - H3N2",
                  "Influenza A - H1N1", "Influenza B - sin linaje", "Influenza B - linaje Victoria")

clasif_controles <- c("Estudiados por pruebas moleculares sin resultados positivos",
                      "Virus sincicial respiratorio (VSR)","SARS-COV-2 por métodos moleculares")
meses <- c("2022_5", "2022_6",  "2022_7",  "2022_8",  "2022_9", "2022_10", "2022_11")

eliminar_sexo_no_binario = T
eliminar_edades_nulas = T




### LECTURA Y PRE-PROCESAMIENTO COMÚN #####
setwd("/home/santiago/Documents/MinSalNacion/influenza/bases_segundo_informe/")
vac_path <- "SNVS_VR_INTERNADO_FALLECIDO_VACUNAS.csv"
data <- read.csv(vac_path)
# Esto mete duplicados. Importante, no debería salvo con refuerzo.

vac_data <- data %>% dplyr::select(CODIGO_CIUDADANO,
                            F_INFLUENZA,
                            F_PRIMERA,
                            F_SEGUNDA,
                            F_ADICIONAL) %>% distinct()

inf_data <- read.csv("Base_Respi_Internado - Base_Respi_Internado.csv")

# Me quedo con los casos que fueron estudiados.
flu_data <- inf_data %>% filter(ESTUDIADO_INFLUENZA == "Si")

# Junto con la información de vacunación de la otra base.
flu_data <- left_join(flu_data, vac_data, by = "CODIGO_CIUDADANO")

# Agrego un campo de fecha mínima para trabajar con una sola
flu_data <- flu_data%>% mutate(FECHA_EST = coalesce(FIS,FECHA_CONSULTA,FTM, FECHA_INTERNACION,FECHA_APERTURA))

# Cambio el formato de fechas para trabajar fácil.
flu_data$FECHA_ESTUDIO <- as.numeric(as.Date(flu_data$FECHA_EST,tryFormats = c("%m/%d/%Y")), origin = "1970-01-01")

flu_data$F_INFLUENZA <-  as.numeric(as.Date(as.POSIXct(flu_data$F_INFLUENZA, format="%m/%d/%Y  %H:%M")))
flu_data$F_PRIMERA <-  as.numeric(as.Date(as.POSIXct(flu_data$F_PRIMERA, format="%m/%d/%Y  %H:%M")))
flu_data$F_SEGUNDA <-  as.numeric(as.Date(as.POSIXct(flu_data$F_SEGUNDA, format="%m/%d/%Y  %H:%M")))
flu_data$F_ADICIONAL <-  as.numeric(as.Date(as.POSIXct(flu_data$F_ADICIONAL, format="%m/%d/%Y  %H:%M")))

# Defino estados de vacunación según diferencia entre la fecha de estudio y las vacunaciones.
flu_data <- flu_data %>% mutate(
  fluvac = case_when(F_INFLUENZA > as.numeric(as.Date("2022-01-01")) & F_INFLUENZA+14 < FECHA_ESTUDIO~'Vac2022',
                                      F_INFLUENZA > as.numeric(as.Date("2022-01-01")) & F_INFLUENZA < FECHA_ESTUDIO~'Vac2022_14dias_previos',
                                      F_INFLUENZA < FECHA_ESTUDIO~"VacPre2022", TRUE ~ "NoVac"),
  covidvac_pre_estudio = case_when(F_PRIMERA+14 < FECHA_ESTUDIO & F_SEGUNDA+14 < FECHA_ESTUDIO & F_ADICIONAL+14 < FECHA_ESTUDIO~3,
                                   (F_PRIMERA+14 < FECHA_ESTUDIO & F_SEGUNDA+14 < FECHA_ESTUDIO)~2,
                                   F_PRIMERA+14 < FECHA_ESTUDIO~1,
                                   TRUE ~0 ),
  MES_ESTUDIO = month(as.Date(FECHA_ESTUDIO, origin = "1970-01-01")),
  AÑO_ESTUDIO = year(as.Date(FECHA_ESTUDIO, origin = "1970-01-01")),
  MES_AÑO_ESTUDIO = paste(AÑO_ESTUDIO,MES_ESTUDIO, sep = "_")
  
  )
# Elimino a los que se vacunaron justo antes de enfermar
flu_data <- flu_data %>% filter((fluvac != "Vac2022_14dias_previos"))

# Detección y eliminación de registros duplicados
flu_data <- flu_data %>%
  group_by(CODIGO_CIUDADANO) %>%
  mutate(n_reg = mean(row_number())) %>% filter(n_reg == 1)

################## CRITERIOS DE SELECCIÓN DE CASOS Y CONTROLES ###############
# Me quedo con los influenza
casos <- flu_data %>% filter(CLASIFICACION_ALGORITMO %in% clasif_casos)
# Elimino codetecciones

casos$status <- "caso"

# Me quedo con los negativo en todo como primer grupo control
controles <- flu_data %>% filter(CLASIFICACION_ALGORITMO %in% clasif_controles)

controles$status <- "control"

case_control <- rbind(casos, controles)
case_control$status <- ifelse(case_control$status == "caso",1,0)

################################# FILTRADO #####################################
# Me quedo con los casos con sexo registrado, aunque se puede corregir manual.
if (eliminar_sexo_no_binario){
  case_control <- case_control %>% filter(SEXO %in% c("M","F"))  
}
# Elimino casos sin edad.
if (eliminar_edades_nulas){
  case_control <- case_control %>% filter(!is.na(EDAD_DIAGNOSTICO))  
}

# Voy a eliminar casos de antes de abril, porque son muy pocos.
# Ojo, cuando aparezcan meses de dos cifras esto se va a romper.
case_control <- case_control %>% filter(MES_AÑO_ESTUDIO %in% meses)

file_name <- paste("flu_ve",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),
                   ".csv",sep = "_")
case_control %>% select(EDAD_DIAGNOSTICO,SEXO, PROVINCIA_RESIDENCIA, covidvac_pre_estudio,fluvac,
                        FALLECIDO, FECHA_ESTUDIO, MES_ESTUDIO, MES_AÑO_ESTUDIO,EVENTO, CLASIFICACION_ALGORITMO, status) %>% 
              write.csv(file_name)


