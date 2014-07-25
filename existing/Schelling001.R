# 24/9/2013
# Jon Minton
# Attempt to create a simple Schelling model in R

rm(list=ls())


#######################################################################################################
# FUNCTIONS

CalcDissimilarity <- function(X, minval, majval, aunits=10){
  # X is the matrix about which the dissimilarity should be calculated
  # minval: the value indicating the presence of a minority agent
  # majval: the value indicating the presence of a majority agent
  
  yband <- dim(X)[1] / aunits
  xband <- dim(X)[2] / aunits
  
  # global measure
  B <- length(which(as.vector(X)==minval)) / length(as.vector(X))
  W <- length(which(as.vector(X)==majval)) / length(as.vector(X))
  
  # local measure
  b <- vector("numeric", aunits)
  w <- vector("numeric", aunits)
  
  for (i in 1:aunits){
  
    x <- X[1:yband + (i - 1) * yband, 1:xband + (i - 1) * xband]
    b[i] <- length(which(as.vector(x)==minval))/length(as.vector(x))
    w[i] <- length(which(as.vector(x)==majval))/length(as.vector(x))  
  }
  
  output <- (1/2) * sum(abs((b/B) - (w/W)))
  
  return(output)
}

InitialiseCells <- function(Nx, Ny, p.red, p.blue){
  output <- matrix(0, nrow=Nx, ncol=Ny)
  
  positions.allowed <- Nx * Ny
  
  n.red   <- as.integer(positions.allowed * p.red)
  n.blue  <- as.integer(positions.allowed * p.blue)
  
  posval.red <- sample((0:(positions.allowed-1)), n.red)
  # can't be in the same positions
  posval.blue <- sample((0:(positions.allowed-1))[!((0:(positions.allowed-1)) %in% posval.red)], n.blue)
  
  xpos.red <- (posval.red %% Nx) + 1
  ypos.red <- (posval.red %/% Ny) + 1
  
  xpos.blue <- (posval.blue %% Nx) + 1 
  ypos.blue <- (posval.blue %/% Ny) + 1 
  
  pos.red <- list(x=xpos.red, y=ypos.red)
  pos.blue <- list(x=xpos.blue, y=ypos.blue)
  
  
  # add reds (value 1)
  
  for (i in 1:n.red){
    output[pos.red$x[i], pos.red$y[i]] <- 1
  }
  
  for (i in 1:n.blue){
    output[pos.blue$x[i], pos.blue$y[i]] <- 2
  }
  
  return(output)
}

CalcNUnlikeCells <- function(X){
  this.celltype <- X[2,2]
  neighbours <- X[-c(2,2)]
  neighbours <- neighbours[neighbours!=0] # remove empty cells
  n.others <- length(which(neighbours!=this.celltype))
  return(n.others)
}

FindEmptyCells <- function(X){
# X is the neighbourhood
  output <- c() # empty vector, populated as new options identified
  
  if (X[1,1]==0){ output <- c(output, "NW")}
  if (X[2,1]==0){ output <- c(output, "W")}
  if (X[3,1]==0){ output <- c(output, "SW")}
  if (X[1,2]==0){ output <- c(output, "N")}
  if (X[3,2]==0){ output <- c(output, "S")}
  if (X[1,3]==0){ output <- c(output, "NE")}
  if (X[2,3]==0){ output <- c(output, "E")}
  if (X[3,3]==0){ output <- c(output, "SE")}
  
  return(output)
}

MakeMove <- function(input){
  # input is output from FindEmptyCells
  # choose one of the possible places to move
  choice <- input[sample(1:length(input), 1)]
  
  if (choice=="NW"){ xadj <- -1 ; yadj <- -1 } 
  if (choice=="W") { xadj <- -1 ; yadj <-  0 }
  if (choice=="SW"){ xadj <- -1 ; yadj <-  1 }
  if (choice=="N") { xadj <-  0 ; yadj <- -1 }
  if (choice=="S") { xadj <-  0 ; yadj <-  1 }
  if (choice=="NE"){ xadj <-  1 ; yadj <- -1 }    
  if (choice=="E") { xadj <-  1 ; yadj <-  0 }
  if (choice=="SE"){ xadj <-  1 ; yadj <-  1 }
  
  output <- list(choice=choice, x=xadj, y=yadj)
  return(output)
}

# ChooseMovetoCell <- function(x){
#   # x is a vector giving which adjacent cells can be moved to
#   # 1 4 6
#   # 2   7
#   # 3 5 8
#   
#   next.pos <- sample(x, 1) # pick one of the possible cells at random
# 
#   if (next.pos==1){
#     i.adj <- -1
#     j.adj <- -1
#   } else if (next.pos==2){
#     i.adj <- -1
#     j.adj <- 0
#   } else if (next.pos==3){
#     i.adj <- 1
#     j.adj <- -1
#   } else if (next.pos==4){
#     i.adj <- -1
#     j.adj <- 0
#   } else if (next.pos==5){
#     i.adj <- 1
#     j.adj <- 1    
#   } else if (next.pos==6){
#     i.adj <- -1
#     j.adj <- 1
#     
#   } else if (next.pos==7){
#     i.adj <- 0
#     j.adj <- 1
#     
#   } else if (next.pos==8){
#     i.adj <- 1
#     j.adj <- 1
#     
#   }
#   output <- list(i.adj=i.adj, j.adj=j.adj)
#   return(output)  
# }

ShuffleCells <- function(inputmatrix, tau, debug=F){
  outputmatrix <- inputmatrix
  Ny <- dim(inputmatrix)[1]
  Nx <- dim(inputmatrix)[2]
  
  for (i in 2:(Ny-1)){
    for (j in 2:(Nx-1)){
      # The neighbourhood comprises the cell of interest ('me')
      # and the eight neighbouring cells
      neighbourhood <- outputmatrix[(i-1):(i+1), (j-1):(j+1)]
      # 'me' is the centre of the neighbourhood
      me <- neighbourhood[2,2]
      # 'Neighbours' are the eight surrounding cells
      # removing the central cell from the 3x3 matrix turns
      # the output into a vector
      neighbours <- neighbourhood[-c(2,2)]
      # the value '1' refers to a red cell, and 
      # '2' to a blue cell
      # The code below counts how many of neighbours are either
      # red or blue
      n.red.neighbours <- length(which(neighbours==1))
      n.blue.neighbours <- length(which(neighbours==2))
      n.empty.neighbours <- length(which(neighbours==0))
      
      if (me==1){ # if I am red
        if (n.blue.neighbours > tau){
          # desire to move
          if (debug==T){
            print(outputmatrix)
            print(c(i, j))
            print(neighbourhood)
            print("I am red and want to move")
            browser()
          }
          # Need to check whether the agent can move
          if (n.empty.neighbours > 0){
            moveoptions <- FindEmptyCells(neighbourhood)
            if (debug==T){
              print("I can move")
              print(moveoptions)
              browser()              
            }
            #move to one of the possible options
            
            mve <- MakeMove(moveoptions)            
            # empty the central cell in the neighbourhood 
            outputmatrix[i,j] <- 0
            # move 'me' to the chosen neighbouring cell

            if (outputmatrix[i+mve$y, j+mve$x]==0){
              outputmatrix[i+mve$y, j + mve$x] <- me
            } else {
              break("Trying to move into occupied cell")
            }

            
            if (debug==T){
              print("I have moved")
              print(mve$choice)
              print(outputmatrix)
              browser()              
            }    
          } else {
            if (debug==T){
              print("I cannot move")
              browser()
              
            }
            
          }
          # 
        }
        
      }
      
      if (me==2){ # if I am blue
        if (n.red.neighbours > tau){
          # desire to move
          if (debug==T){
            print(outputmatrix)
            print(c(i, j))
            print(neighbourhood)
            print("I am blue and want to move")
            browser()            
          }
          if (n.empty.neighbours > 0){
            # I can move
            moveoptions <- FindEmptyCells(neighbourhood)
            if (debug==T){
              print("I can move")
              print(moveoptions)
              browser()
            } 
            #move to one of the possible options
            
            mve <- MakeMove(moveoptions)            
            # empty the central cell in the neighbourhood 
            outputmatrix[i,j] <- 0
            # move 'me' to the chosen neighbouring cell
            outputmatrix[i+mve$y, j+mve$x] <- me
            if (debug==T){
              print("I have moved")
              print(mve$choice)
              print(outputmatrix)
              browser()              
            }
            
            
            
          } else {
            if (debug==T){
              print("I cannot move")
              browser()
              
            }
            
          }
          
        }
        
      }
    }
    
  }
  
  return(outputmatrix)
}


    
#     
#     # Only move if 
#     #  1) central cell occupied
#     #  2) At least one neighbouring unoccupied cell to move to
#     #  3) the number of unlike neighbours exceeds the threshold
#     # Otherwise, stay in position
#     
#     if (this.cell != 0){
#       # stay or move, depending on n unlike cells
#       n.unlikecells <- CalcNUnlikeCells(localcells)
#       
#       if (n.unlikecells >= tau){
#         # the agent wants to move, is there somewhere to move to?
#         emptycells <- CalcNSpareCells(localcells)
#         n.emptycells <- length(emptycells)
#         if(n.emptycells!=0){
#           # a desire to move, and opportunity to move: which cell is moved to?
#           m <- ChooseMovetoCell(emptycells)
#           dtablock[i - m$i.adj, j - m$j.adj, k] <- this.cell
#           dtablock[i,j,k] <- 0 # moved from
#         }
#       }  
#     }
#     
#     if (debug==T){
#       image(dtablock[,,k], xaxt="n", yaxt="n", col=c("white","red","blue"))
#       print(paste(i, " ", j, " ", k))
#       browser()
#     }
#     
#   }
# }

#########################################################################################################
# INPUTS 


# Input params

# matrix dimensions
# proportion reds 
# proportion blues
# iterations 
dbg <- T

p.red <- 0.35
p.blue <- 0.35

Nx <- Ny <- 100 # width and height
K <- 200 # iterations

tau <- 2 # move if four or more neighbours are of other type

cellmatrix <- InitialiseCells(Nx, Ny, p.red, p.blue)
# creates a block 50 rows tall, 50 cols wide, and 50 deep

image(t(cellmatrix), xaxt="n", yaxt="n", col=c("white", "red", "blue"))

# now to 'move' agents depending on value of adjacent cells

# to begin with, start from 2nd col and row, going to i-1th col and row, to avoid
# having to think about edge effects for now

# initially, cells on boundaries don't move (improve on this later)

cellmatrix.list <- vector("list", length=K)

# k: number of iterations
# i : width
# j : height
cellmatrix.list[[1]] <- cellmatrix
for (k in 2:K){
#  print(cellmatrix)
  cellmatrix.old <- cellmatrix
  cellmatrix <- ShuffleCells(cellmatrix, tau, debug=F) 
  cellmatrix.list[[k]] <- cellmatrix
  cellmatrix.dif <- cellmatrix - cellmatrix.old
}
       

# Check no agents are disappearing 

sapply(cellmatrix.list, function(x) sum(as.vector(x)))



for (i in 1:K){
  image(t(cellmatrix.list[[i]]), xaxt="n", yaxt="n", main=i, col=c("white", "red", "blue"))
  browser()
}


# to do: 
# Adapt to two state (correct) version
# Remove boundary issues by turning area into a torus




#Sys.sleep(0.1)
# 