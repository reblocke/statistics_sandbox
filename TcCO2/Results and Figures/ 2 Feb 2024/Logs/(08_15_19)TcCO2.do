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

