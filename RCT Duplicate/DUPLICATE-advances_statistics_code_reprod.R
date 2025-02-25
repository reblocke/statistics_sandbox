#
# -- Advances_Statistics_Code_Reprod.R: statistical facets of reproducibility --
#
rm      ( list = ls ( ) )                                  #  remove all data objects

library ( beeswarm )
library ( ggplot2  )


# -- Define defaults for graphics files ----------------------------------------
#
Save.Graphics       <- FALSE

#    Save.Graphics <- TRUE                                       #  if TRUE, script saves graphics to ps files
#
ps.options ( onefile   = FALSE, paper = 'special', height = 7, width = 7 )
options    ( papersize = 'special' )


# -- Define population parameters, alpha, and sample numbers -------------------
#
PopMean      <-    0                                       #  population mean
PopSD        <-    1                                       #  population standard deviation
dMean        <-    0.5                                     #  shift in population mean for population 2
nSim         <- 1000                                       #  number of simulations for each sample size

Actual.Alpha <-    0.05

n.1          <-   10                                       #  theoretical power: 0.18
n.2          <-   30                                       #                     0.48
n.3          <-   64                                       #                     0.80
n.Last       <-  100                                       #                     0.94


# -- Create shift for log scale of P values ------------------------------------
#
log.Shift <- 0.0001
Log.Alpha <- log10 ( Actual.Alpha + log.Shift )


# -- Run the simulations -------------------------------------------------------

n.Sample.Sizes <- 4

t.Power        <- matrix ( nrow = n.Sample.Sizes, ncol = 2 )
means.Power    <- matrix ( nrow = n.Sample.Sizes, ncol = 3 )


# -- When the Null Hypothesis is False :: the simulation -----------------------                      # -- the simulation: first line
#
if ( Save.Graphics == TRUE ) postscript ( file = 'Sample_Observations.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1 )

for ( n in 1:n.Sample.Sizes )
    { if ( n == 1              ) nObs <- n.1
      if ( n == 2              ) nObs <- n.2
      if ( n == 3              ) nObs <- n.3
      if ( n == n.Sample.Sizes ) nObs <- n.Last

      Simulation.Values <- matrix ( nrow = nSim, ncol = 3 )

      t.Power     [ n, 1 ] <- nObs
      means.Power [ n, 1 ] <- nObs

      for ( s in 1:nSim )
          { The_Data  <- matrix ( nrow = 2 * nObs, ncol = 2 )

            Simulation.Number          <- s
            Simulation.Values [ s, 1 ] <- Simulation.Number

            for ( k in 1:nObs )                                                                       #  group 0 obs
                { The_Data [ k,1 ] <- 0
                  The_Data [ k,2 ] <- round ( rnorm ( 1, mean = PopMean, sd = PopSD ), 3 )
                  }

            for ( m in ( 1 + nObs ) : ( nObs + nObs ) )                                               #  group 1 obs
                { The_Data [ m,1 ] <- 1
                  The_Data [ m,2 ] <- round ( rnorm ( 1, mean = PopMean + dMean, sd = PopSD ), 3 )
                  }

            tTest.P <- t.test ( The_Data [ ,2 ] ~ The_Data [ ,1 ],                                    #  2-sample t test
                                alternative = c ( 'two.sided' ),
                                var.equal   = FALSE,
                                conf.level  = 1 - Actual.Alpha
                                )

            Simulation.Values [ s,2 ] <- round ( tTest.P$p.value, 4 )
            Simulation.Values [ s,3 ] <- abs ( tTest.P$estimate [ 2 ] - tTest.P$estimate [ 1 ] )      #  estimated diff in means
            d.Mu                      <- round ( Simulation.Values [ s,3 ], 3 )


            if ( n == 1 && s < 4 )
               { beeswarm ( The_Data [ ,2 ] ~ The_Data [ ,1 ],
                            main = paste ( 'Sim', s, ':: dmu =', d.Mu, ':: P =', round ( tTest.P$p.value, 3 ) ),
                            ylab = 'Diff in means',
                            ylim = c ( -3, +3 ),
                            xlab = '',
                            pch  = 19,
                            col  = c ( 'gray50', 'black' ),
                            bty  = 'n'
                            )
                 #
                 lines ( x = c ( 0.85, 1.15 ), y = rep ( mean ( The_Data [ ,2 ] [ The_Data [ ,1 ] == 0 ] ), 2 ) )
                 lines ( x = c ( 1.85, 2.15 ), y = rep ( mean ( The_Data [ ,2 ] [ The_Data [ ,1 ] == 1 ] ), 2 ) )
                 #
                 if ( Save.Graphics == FALSE ) Sample_Observations <- recordPlot ( )
                 #
                 }

            }    #  for ( s in 1:nSim )


      # -- Store results from each sample size run -----------------------------
      #
      if ( n == 1              ) Simulation.Values.Sample.Size.1    <- Simulation.Values
      if ( n == 2              ) Simulation.Values.Sample.Size.2    <- Simulation.Values
      if ( n == 3              ) Simulation.Values.Sample.Size.3    <- Simulation.Values
      if ( n == n.Sample.Sizes ) Simulation.Values.Sample.Size.Last <- Simulation.Values

      if ( n == 1 ) dMax <- max ( Simulation.Values [ ,3 ] )

      t.Power     [ n, 2 ] <- mean ( Simulation.Values [ ,2 ] < Actual.Alpha )
      means.Power [ n, 2 ] <- mean ( Simulation.Values [ ,3 ] > dMean )


      # -- Get proportion of delta means > dMean when P < alpha ----------------
      #
      means.Power [ n, 3 ] <- round ( mean ( Simulation.Values [ ,3 ] [ Simulation.Values [ ,2 ] < Actual.Alpha ] > dMean ), 3 )

      }    #  for ( n in 1:n.Sample.Sizes )

if ( Save.Graphics == TRUE ) dev.off ( )


# -- Print power for P value and delta mean > dMean ----------------------------
#
t.Power        #  [ ,1 ]: sample size each of 2 groups
               #  [ ,2 ]: proportion of simulations in which P < alpha

means.Power    #  [ ,1 ]: sample size each of 2 groups
               #  [ ,2 ]: proportion of simulations in which delta mean > dMean
               #  [ ,3 ]: proportion of simulations in which delta mean > dMean for those simulations in which P < alpha


# -- Generate distribution of log P --------------------------------------------
#
P.Value.Percentiles <- matrix ( nrow = n.Sample.Sizes, ncol = 5 )


if ( Save.Graphics == TRUE ) postscript ( file = 'Log10_P_Histograms.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1 )

for ( n in 1:n.Sample.Sizes )
    { if ( n == 1              ) { Simulation.Values <- Simulation.Values.Sample.Size.1
                                   nObs              <- n.1
                                   }
      if ( n == 2              ) { Simulation.Values <- Simulation.Values.Sample.Size.2
                                   nObs              <- n.2
                                   }
      if ( n == 3              ) { Simulation.Values <- Simulation.Values.Sample.Size.3
                                   nObs              <- n.3
                                   }
      if ( n == n.Sample.Sizes ) { Simulation.Values <- Simulation.Values.Sample.Size.Last
                                   nObs              <- n.Last
                                   }

      P.Value.Percentiles [ n, 1 ] <- n
      P.Value.Percentiles [ n, 2 ] <- nObs
      P.Value.Percentiles [ n, 3 ] <- mean ( Simulation.Values [ ,2 ] < 0.001        ) * nSim
      P.Value.Percentiles [ n, 4 ] <- mean ( Simulation.Values [ ,2 ] < 0.01         ) * nSim
      P.Value.Percentiles [ n, 5 ] <- mean ( Simulation.Values [ ,2 ] < Actual.Alpha ) * nSim

      Simulation.Values [ ,2 ] <- log10 ( Simulation.Values [ ,2 ] + log.Shift )

      hist ( Simulation.Values [ ,2 ],
             freq   = FALSE,
             main   = paste ( 'n =', nObs, ':: log P values' ),
             xlim   = c ( log10 ( log.Shift ), log10 ( 1 + log.Shift ) ),
             xlab   = 'log10 P value',
             ylab   = '',
             ylim   = c ( 0, 1.5 ),
             nclass = 20    #  20 .. 24 .. 20 best
             )
      lines  ( density ( Simulation.Values [ ,2 ], adjust = 1, na.rm = TRUE ), col = 'blue', lwd = 1 )
      segments ( Log.Alpha, 0, Log.Alpha, 10, lwd = 1, col = 'red' )

      }    #  for ( n in 1:n.Sample.Sizes )

if ( Save.Graphics == TRUE ) dev.off ( )
if ( Save.Graphics == FALSE ) Log.P.Histograms <- recordPlot ( )


# -- Print out numbers of P values less than cutoffs of 0.001, 0.01, and 0.05 --
#
colnames ( P.Value.Percentiles ) <- c ( 'Sample', 'nObs', 'n_0.001', 'n_0.01', 'n_0.05' )
#
P.Value.Percentiles


# -- Generate distribution of estimated difference in means --------------------
#
dMax <- ceiling ( dMax )                                                                              #  maximum estimated difference in means

if ( Save.Graphics == TRUE ) postscript ( file = 'Delta_Mu_Histograms.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1 )

for ( n in 1:n.Sample.Sizes )
    { if ( n == 1              ) { Simulation.Values <- Simulation.Values.Sample.Size.1
                                   nObs              <- n.1
                                   }
      if ( n == 2              ) { Simulation.Values <- Simulation.Values.Sample.Size.2
                                   nObs              <- n.2
                                   }
      if ( n == 3              ) { Simulation.Values <- Simulation.Values.Sample.Size.3
                                   nObs              <- n.3
                                   }
      if ( n == n.Sample.Sizes ) { Simulation.Values <- Simulation.Values.Sample.Size.Last
                                   nObs              <- n.Last
                                   }

      Temp.dMean.Greater <- subset ( Simulation.Values, Simulation.Values [ ,2 ] < Actual.Alpha )


      # -- Distribution of delta mu --------------------------------------------
      #
      hist ( Temp.dMean.Greater [ ,3 ],
             freq   =  FALSE,
             main   =  paste ( 'n =', nObs, ':: delta mu' ),
             xlim   = c ( 0, dMax ),
             xlab   = 'Delta mu',
             ylab   = '',
             ylim   = c ( 0, 3.5 ),
             nclass = 20
             )
      lines  ( density ( Temp.dMean.Greater [ ,3 ], adjust = 1, na.rm = TRUE ), col = 'blue', lwd = 1 )
      abline ( v = dMean, lwd = 1, col = 'red' )                                                      #  true difference in means

      }    #  for ( n in 1:n.Sample.Sizes )
#
if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Delta.Mu.Histograms <- recordPlot ( )


# -- GGPlot density plots ------------------------------------------------------
#
if ( Save.Graphics == TRUE ) postscript ( file = 'Contour_Plot_Sample_Size_%03d.ps' )

for ( n in 1:n.Sample.Sizes )
    { if ( n == 1              ) { Simulation.Values.3D.Density <- Simulation.Values.Sample.Size.1
                                   nObs                         <- n.1
                                   }
      if ( n == 2              ) { Simulation.Values.3D.Density <- Simulation.Values.Sample.Size.2
                                   nObs                         <- n.2
                                   }
      if ( n == 3              ) { Simulation.Values.3D.Density <- Simulation.Values.Sample.Size.3
                                   nObs                         <- n.3
                                   }
      if ( n == n.Sample.Sizes ) { Simulation.Values.3D.Density <- Simulation.Values.Sample.Size.Last
                                   nObs                         <- n.Last
                                   }

      yMin <- round ( min ( Simulation.Values.3D.Density [ ,3 ] [ Simulation.Values.3D.Density [ ,2 ] < Actual.Alpha ] ), 2 )
      yMax <- round ( max ( Simulation.Values.3D.Density [ ,3 ] [ Simulation.Values.3D.Density [ ,2 ] < Actual.Alpha ] ), 2 )


      Density.Graphic <- ggplot ( data.frame ( Simulation.Values.3D.Density ),
                                  aes        ( x = Simulation.Values.3D.Density [ ,2 ] + log.Shift,
                                               y = Simulation.Values.3D.Density [ ,3 ]
                                               )
                                  ) +
                                  theme ( legend.position = 'none' ) +

                                  scale_x_log10   ( limits = c ( log.Shift, 1.0000 ) ) +
                                  coord_cartesian ( ylim = c ( 0, dMax ) ) +
                                  xlab    ( 'log10 P' ) +
                                  ylab    ( '|Diff in means|' ) +
                                  ggtitle ( paste ( n, ' ... n =', nObs, '::', sprintf ( '%.2f', yMin ), 'to', sprintf ( '%.2f', yMax ) ) ) +

                                  stat_density2d      ( aes ( fill = ..level.. ), geom = 'polygon' ) +
                                  scale_fill_gradient ( low = 'gray80', high = 'darkslategray' ) +

                                  geom_point          ( col = 'firebrick4', size = 1.5, alpha = 1.000 ) +

                                  geom_density2d      ( col = 'white' ) +

                                  geom_vline          ( xintercept = 0.001, col = 'steelblue' ) +
                                  geom_vline          ( xintercept = 0.05,  col = 'steelblue' ) +

                                  geom_hline          ( yintercept = dMean, col = 'steelblue' ) +
                                  geom_hline          ( yintercept = yMin,  col = 'purple'    ) +
                                  geom_hline          ( yintercept = yMax,  col = 'purple'    )

      print ( Density.Graphic )

      if ( Save.Graphics == FALSE ) { if ( n == 1              ) All.P.Densities.1    <- recordPlot ( )
                                      if ( n == 2              ) All.P.Densities.2    <- recordPlot ( )
                                      if ( n == 3              ) All.P.Densities.3    <- recordPlot ( )
                                      if ( n == n.Sample.Sizes ) All.P.Densities.Last <- recordPlot ( )
                                      }
      }

if ( Save.Graphics == TRUE ) dev.off ( )
#
# ------------------------------------------------------------------------------                      # -- the simulation: last line


# -- Replay each data graphic --------------------------------------------------
#
if ( Save.Graphics == FALSE ) replayPlot ( Sample_Observations  )    #  Figure 2: first 4 sets of sample observations
if ( Save.Graphics == FALSE ) replayPlot ( Log.P.Histograms     )    #  Figure 3: distributions of log P
if ( Save.Graphics == FALSE ) replayPlot ( Delta.Mu.Histograms  )    #  Figure 3: distributions of estimated difference in means

if ( Save.Graphics == FALSE ) replayPlot ( All.P.Densities.1    )    #  Figure 3: sample size 1: density plot
if ( Save.Graphics == FALSE ) replayPlot ( All.P.Densities.2    )    #  Figure 3: sample size 2: density plot
if ( Save.Graphics == FALSE ) replayPlot ( All.P.Densities.3    )    #  Figure 3: sample size 3: density plot
if ( Save.Graphics == FALSE ) replayPlot ( All.P.Densities.Last )    #  Figure 3: sample size 4: density plot


# -- Table 1 -------------------------------------------------------------------                      # -- Table 1: first line
#
n.P   <- 32
Alpha <-  0.05                                                                                        # -- critical significance level

Replication.P <- matrix ( nrow = n.P, ncol = 2 )

for ( j in 1:n.P )
    { if ( j ==  1 ) P_1 <- .0122274 #leader
      if ( j ==  2 ) P_1 <- .0055911 #declare timi
      if ( j ==  3 ) P_1 <- .0418784 #empa-reg
      if ( j ==  4 ) P_1 <- .0213696 #canvas
      if ( j ==  5 ) P_1 <- .789009 #carmalina
      if ( j ==  6 ) P_1 <- .7244875 #tecos
      if ( j ==  7 ) P_1 <- 1 #SAVOR TIMI
      if ( j ==  8 ) P_1 <- 1 #LEAD 2
      if ( j ==  9 ) P_1 <- .0000909 #Triton TIMI
      if ( j == 10 ) P_1 <- .0001381 #PLATO
      if ( j == 11 ) P_1 <- .0067104 #ISAR-REACT
      if ( j == 12 ) P_1 <- .0111473 #ARIsTOTLE
      if ( j == 13 ) P_1 <- .0002096 #RELY
      if ( j == 14 ) P_1 <- .0135936 #ROCKET-AF
      if ( j == 15 ) P_1 <- .0784689 #EISTEIN DVT
      if ( j == 16 ) P_1 <- .5938132 #EISTEIN PE
      if ( j == 17 ) P_1 <- .7830252 #RECOVER 1
      if ( j == 18 ) P_1 <- .3168379 #AMPLIFY
      if ( j == 19 ) P_1 <- 9.23e-06 #RECORDI
      if ( j == 20 ) P_1 <- .2095088 #TRANSCEND
      if ( j == 21 ) P_1 <- .8043096 #ONTARGET
      if ( j == 22 ) P_1 <- .0024478 #HORIZON PFT
      if ( j == 23 ) P_1 <- .0001769 #VERO
      if ( j == 24 ) P_1 <- 3.52e-08 #DAPA CKD
      if ( j == 25 ) P_1 <- 9.05e-07 #PARADIGM HF
      if ( j == 26 ) P_1 <- 5.18e-06 #P04334Af
      if ( j == 27 ) P_1 <- .7698346 #D5896
      if ( j == 28 ) P_1 <- 1.07e-07 #1M PACTe'g
      if ( j == 29 ) P_1 <- 3.82e-06 #POET-COPD
      if ( j == 30 ) P_1 <- .691256 #INSPIRE
      if ( j == 31 ) P_1 <- .8074158 #CAROLINA
      if ( j == 32 ) P_1 <- .5444627 #Pronounce
      
      z.Crit <- qnorm ( 1 - ( Alpha / 2 ) )                                                           # -- critical value of z associated with Alpha
      z_1    <- qnorm ( 1 - ( P_1 / 2 ) )                                                             # -- observed value of z associated with initial value of P

      Pr.P.Duplicate <- pnorm ( z.Crit - z_1, lower.tail = FALSE )                                    # -- probability a duplicate experiment will achieve P < Alpha

      Replication.P [ j, 1 ] <- P_1
      Replication.P [ j, 2 ] <- round ( Pr.P.Duplicate, 5 )
      }
#
colnames ( Replication.P ) <- c ( 'P_1', 'Pr{P_2<0.05}' )
#
Replication.P
#
# ------------------------------------------------------------------------------                      # -- Table 1: last line

