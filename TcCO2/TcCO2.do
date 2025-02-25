clear
cd "/Users/blocke/Box Sync/Residency Personal Files/Stats/TcCO2"
//cd "C:\Users\reblo\Box\Residency Personal Files\Stats\Clinical Trials Processed"


program define datetime 
end

capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "TcCO2.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"

set scheme cleanplots
graph set window fontface "Helvetica"



capture log close
log using "Results and Figures/$S_DATE/Logs/temp.log", append

import excel "TcCo2 data table.xlsx", sheet("Sheet1") firstrow case(lower)


twoway (scatter carbondioxidetensionbefore oxygensaturationbeforestudy) (lfit carbondioxidetensionbefore oxygensaturationbeforestudy), title("Before Treatment") legend(off) 

twoway (scatter carbondioxidetensionafter oxygensaturationafterstudy) (lfit carbondioxidetensionbefore oxygensaturationafterstudy), title("After Treatment") legend(off) 


gen delta_ph = phafterstudy - phbeforestudy
gen delta_co2 = carbondioxidetensionafter - carbondioxidetensionbefore
gen delta_po2 = oxygentensionafterstudy - oxygentensionbeforestudy
gen delta_base = baseexcessafterstudy - baseexcessbeforestudy


twoway (scatter delta_co2 delta_base) (lfit delta_co2 delta_base), title("Delta") legend(off) 
twoway (scatter delta_co2 delta_po2) (lfit delta_co2 delta_po2), title("Delta") legend(off) 



/* ANCOVA ? */ 

regress carbondioxidetensionafter carbondioxidetensionbefore ageyrs bmi
