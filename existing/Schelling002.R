# 25/9/2013
# Jon Minton
# Attempt to create a simple Schelling model in R
# Need to create the correct, two state, version of the model



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
  B <- length(which(as.vector(X)==minval))
  W <- length(which(as.vector(X)==majval)) 
  
  # local measure
  b <- vector("numeric", aunits)
  w <- vector("numeric", aunits)
  
  for (i in 1:aunits){
    
    x <- X[1:yband + (i - 1) * yband, 1:xband + (i - 1) * xband]
    b[i] <- length(which(as.vector(x)==minval))
    w[i] <- length(which(as.vector(x)==majval))  
  }
  
  output <- (1/2) * sum(abs((b/B) - (w/W)))
  
  return(output)
}

InitialiseCells <- function(Nx, Ny, p.black){
  output <- matrix(0, nrow=Nx, ncol=Ny)
  
  positions.allowed <- Nx * Ny
  
  n.black   <- as.integer(positions.allowed * p.black)
  
  posval.black <- sample((0:(positions.allowed-1)), n.black)
  # can't be in the same positions
  
  xpos.black <- (posval.black %% Nx) + 1
  ypos.black <- (posval.black %/% Ny) + 1
  
  
  pos.black <- list(x=xpos.black, y=ypos.black)
  
  # add reds (value 1)
  
  for (i in 1:n.black){
    output[pos.black$x[i], pos.black$y[i]] <- 1
  }
    
  return(output)
}

CalcNUnlikeCells <- function(X){
  this.celltype <- X[2,2]
  neighbours <- X[-c(2,2)]
  n.others <- length(which(neighbours!=this.celltype))
  return(n.others)
}

FindSwapCells <- function(X, mycol){
# X is the neighbourhood
  output <- c() # empty vector, populated as new options identified
  
  if (X[1,1]!=mycol){ output <- c(output, "NW")}
  if (X[2,1]!=mycol){ output <- c(output, "W")}
  if (X[3,1]!=mycol){ output <- c(output, "SW")}
  if (X[1,2]!=mycol){ output <- c(output, "N")}
  if (X[3,2]!=mycol){ output <- c(output, "S")}
  if (X[1,3]!=mycol){ output <- c(output, "NE")}
  if (X[2,3]!=mycol){ output <- c(output, "E")}
  if (X[3,3]!=mycol){ output <- c(output, "SE")}
  
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


SwapCells <- function(inputmatrix, tau, debug=F){
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
      # the value '1' refers to a black cell
      # So the value '0' is a white cell
      
      # The code below counts how many of neighbours are either
      # black or white
      n.black.neighbours <- length(which(neighbours==1))
      n.white.neighbours <- length(which(neighbours==0))
      
      if (me==1){ # if I am black
        if (n.white.neighbours > tau){
          # desire to move
          if (debug==T){
            print(outputmatrix)
            print(c(i, j))
            print(neighbourhood)
            print("I am black and want to move")
            browser()
          }
          moveoptions <- FindSwapCells(neighbourhood, 1)
            if (debug==T){
              print("I can move")
              print(moveoptions)
              browser()              
            }
            #move to one of the possible options
            
            mve <- MakeMove(moveoptions)            
            # empty the central cell in the neighbourhood 
            outputmatrix[i,j] <- 0
            outputmatrix[i+mve$y, j + mve$x] <- 1
            
            if (debug==T){
              print("I have swapped")
              print(mve$choice)
              print(outputmatrix)
              browser()              
            }    
          } 
          # 
        }
      if (me==0){ # if I am white
        if (n.black.neighbours > tau){
          # desire to move
          if (debug==T){
            print(outputmatrix)
            print(c(i, j))
            print(neighbourhood)
            print("I am white and want to move")
            browser()            
          }
          
          moveoptions <- FindSwapCells(neighbourhood, 0)
          if (debug==T){
            print("I can move")
            print(moveoptions)
            browser()
          } 
          #move to one of the possible options
          
          mve <- MakeMove(moveoptions)            
          outputmatrix[i,j] <- 1
          # move 'me' to the chosen neighbouring cell
          outputmatrix[i+mve$y, j+mve$x] <- 0
          if (debug==T){
            print("I have moved")
            print(mve$choice)
            print(outputmatrix)
            browser()              
          }
          
          
          
          
          
        }
        
      }
        
      }
      

    }
  
  return(outputmatrix)
}


    


#########################################################################################################
# INPUTS 


# Input params

# matrix dimensions
# proportion reds 
# proportion blues
# iterations 
dbg <- T

p.black <- 0.20

Nx <- Ny <- 100 # width and height
K <- 200 # iterations

tau <- 4 # move if four or more neighbours are of other type

cellmatrix <- InitialiseCells(Nx, Ny, p.black)
# creates a block 50 rows tall, 50 cols wide, and 50 deep

image(t(cellmatrix), xaxt="n", yaxt="n", col=c("white", "black"))

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
#  cellmatrix.old <- cellmatrix
  cellmatrix <- SwapCells(cellmatrix, tau, debug=F) 
#  browser()
  
  cellmatrix.list[[k]] <- cellmatrix
#  cellmatrix.dif <- cellmatrix - cellmatrix.old
}
       

# Check no agents are disappearing 

sapply(cellmatrix.list, function(x) sum(as.vector(x)))



for (i in 1:K){
  image(t(cellmatrix.list[[i]]), xaxt="n", yaxt="n", main=i, col=c("white", "black"))
  browser()
}

dissimilarities <- sapply(cellmatrix.list, function(x) CalcDissimilarity(x, 0, 1))
plot(1:200, dissimilarities, ylab="Dissimilarity index", xlab="Schelling Round")


Df <- data.frame(x=1:200, y=dissimilarities)
nonlinmodel <- nls (dissimilarities ~ Alpha*(1 - exp(-theta * 1:200)), 
                    start=list(Alpha=0.035, theta=0.5)
                    )

lines(1:200, predict(nonlinmodel, 1:200), col="red", lwd=2)

nonlinmodel2 <- nls(dissimilarities ~ c / ( 1 + a * exp(-(b * 1:200))), start=list(a=0.5, b=0.5, c=0.5))

lines(1:200, predict(nonlinmodel2, 1:200), col="blue", lwd=2, lty="dashed")

legend("bottomright", lwd=2, col=c("red", "blue"), lty=c("solid", "dashed"), 
       legend=c(paste("Model 1: AIC=", round(AIC(nonlinmodel), 0)), paste("Model 2: AIC=", round(AIC(nonlinmodel2),0)))
)

# to do: 
# Adapt to two state (correct) version
# Remove boundary issues by turning area into a torus


#######################################################################################################
# Proof of concept : iterating 1000 times

# matrix dimensions
# proportion reds 
# proportion blues
# iterations 

dissimilarities.list <- vector("list", 1000)


for (M in 1:1000){
  p.black <- 0.20
  
  Nx <- Ny <- 100 # width and height
  K <- 200 # iterations
  
  tau <- 4 # move if four or more neighbours are of other type
  
  cellmatrix <- InitialiseCells(Nx, Ny, p.black)
  # creates a block 50 rows tall, 50 cols wide, and 50 deep
  
  image(t(cellmatrix), xaxt="n", yaxt="n", col=c("white", "black"))
  
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
    #  cellmatrix.old <- cellmatrix
    cellmatrix <- SwapCells(cellmatrix, tau, debug=F) 
    #  browser()
    
    cellmatrix.list[[k]] <- cellmatrix
    #  cellmatrix.dif <- cellmatrix - cellmatrix.old
  }
  
  
  # Check no agents are disappearing 
  
  sapply(cellmatrix.list, function(x) sum(as.vector(x)))
  
  dissimilarities.list[[M]] <- sapply(cellmatrix.list, function(x) CalcDissimilarity(x, 0, 1))
}

dismat <- matrix(NA, nrow=200, ncol=176)
for (i in 1:176){
  dismat[,i] <- dissimilarities.list[[i]]  
}



#Sys.sleep(0.1)
# 