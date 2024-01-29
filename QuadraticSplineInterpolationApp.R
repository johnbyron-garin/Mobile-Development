# John Byron Garin
# Exercise 10
# Quadratic Spline Interpolation

library(shiny)

gaussJordan <- function(resultMat){
  lengthofColumn = length(resultMat[,1])
  lengthofRow = length(resultMat[1,])
  
  i = 1
  
  while(i <= lengthofRow-1){
    if((i > 1) & (i != (lengthofRow-1))){
      pivotColumn = resultMat[-c(1:i-1),-c(1:i-1)]
      pivotColumn = pivotColumn[,1]
      z = abs(pivotColumn)
      x = order(z, decreasing = TRUE)
      indexofLargest = x[1]
      indextobeswapwith = i
      largestNum = resultMat[indexofLargest+(i-1),]
      tobeswapwith = resultMat[indextobeswapwith,]
      resultMat[indexofLargest+(i-1),] = tobeswapwith
      resultMat[indextobeswapwith,] = largestNum
    }else{
      pivotColumn = resultMat[,1]
      z = abs(pivotColumn)
      x = order(z, decreasing = TRUE)
      indexofLargest = x[1]
      indextobeswapwith = 1
      largestNum = resultMat[indexofLargest,]
      tobeswapwith = resultMat[indextobeswapwith,]
      resultMat[indexofLargest,] = tobeswapwith
      resultMat[indextobeswapwith,] = largestNum
    }
    # algorithm used to perform the steps of gauss-jordan elimination
    pivotElement = resultMat[i,i]
    # if pivotElement is zero, then it will return NA
    if(pivotElement == 0){
      return (NA)
    }else{
      x = 1
      while(x <= lengthofRow){
        resultMat[i,x] = resultMat[i,x] / pivotElement
        x = x + 1
      }
      pivotRow = resultMat[i,]
      
      j = 1
      while(j <= lengthofColumn){
        # updates the rows except the element that could be found when row==column
        if(j != i){
          multiplier = resultMat[j,i]
          tempVector = pivotRow * multiplier
          resultMat[j,] = resultMat[j,] - tempVector
          j = j + 1
        }else{
          j = j + 1
        }
      }
    }
    i = i + 1
  }
  return(resultMat[,lengthofRow])
}

poly.qsi <- function (data, x){
    dataPoints = length(data[[1]])
    intervals = dataPoints - 1
    equationNum = 3*intervals
    
    colNames <- vector()
    
    Acounter = 3
    sub = 2
    while(Acounter <= equationNum-1){
      colNames[Acounter] = paste("a", as.character(sub), sep = "")
      sub = sub + 1
      Acounter = Acounter + 3
    }
    
    Bcounter = 1
    sub = 1
    while(Bcounter <= equationNum-1){
      colNames[Bcounter] = paste("b", as.character(sub), sep = "")
      sub = sub + 1
      Bcounter = Bcounter + 3
    }
    
    Ccounter = 2
    sub = 1
    while(Ccounter <= equationNum-1){
      colNames[Ccounter] = paste("c", as.character(sub), sep = "")
      sub = sub + 1
      Ccounter = Ccounter + 3
    }
    
    augCoefMatrix <- matrix(
      0,
      nrow = equationNum-1,
      ncol = equationNum-1,
      byrow = TRUE,
      dimnames = list(1:(equationNum-1), colNames)
    )
    
    RHS <- matrix(
      0,
      nrow = equationNum-1,
      ncol = 1,
      byrow = TRUE,
      dimnames = list(1:(equationNum-1), c("RHS"))
    )
    
    equationCounter = 1
    system <- list()
    eq = 1
    
    # condition 1
    i = 3 # 1 indexing
    while(i <= (dataPoints)){
      # equation format 1
      sub = i-1-1
      varA = paste("a", as.character(sub), sep = "")
      coefA = data[[1]][i-1]^2
      termA = paste(as.character(coefA), "*", varA, sep = " ")
      
      if(varA %in% colNames){
        augCoefMatrix[eq, varA] = coefA
      }
      
      sub = i-1-1
      varB = paste("b", as.character(sub), sep = "")
      coefB = data[[1]][i-1]
      termB = paste(as.character(coefB), "*", varB, sep = " ")
      
      if(varB %in% colNames){
        augCoefMatrix[eq, varB] = coefB
      }
      
      sub = i-1-1
      varC = paste("c", as.character(sub), sep = "")
      coefC = 1
      termC = paste(as.character(coefC), "*", varC, sep = " ")
      
      if(varC %in% colNames){
        augCoefMatrix[eq, varC] = coefC
      }
      
      sol = data[[2]][i-1]
      RHS[eq, 1] = sol
      
      
      equation = paste(termA, "+", termB, "+", termC, "=", sol)
      system[[equationCounter]] = equation
      equationCounter = equationCounter + 1
      
      eq = eq + 1
      
      # equation format 2
      sub = i-1
      varA = paste("a", as.character(sub), sep = "")
      coefA = data[[1]][i-1]^2
      termA = paste(as.character(coefA), "*", varA, sep = " ")
      
      if(varA %in% colNames){
        augCoefMatrix[eq, varA] = coefA
      }
      
      sub = i-1
      varB = paste("b", as.character(sub), sep = "")
      coefB = data[[1]][i-1]
      termB = paste(as.character(coefB), "*", varB, sep = " ")
      
      if(varB %in% colNames){
        augCoefMatrix[eq, varB] = coefB
      }
      
      sub = i-1
      varC = paste("c", as.character(sub), sep = "")
      coefC = 1
      termC = paste(as.character(coefC), "*", varC, sep = " ")
      
      if(varC %in% colNames){
        augCoefMatrix[eq, varC] = coefC
      }
      
      sol = data[[2]][i-1]
      RHS[eq, 1] = sol
      
      equation = paste(termA, "+", termB, "+", termC, "=", sol)
      system[[equationCounter]] = equation
      equationCounter = equationCounter + 1
      
      eq = eq + 1
      
      i = i + 1
    }
    
    # condition 2
    
    # equation format 1
    sub = 1
    varA = paste("a", as.character(sub), sep = "")
    coefA = data[[1]][1]^2
    termA = paste(as.character(coefA), "*", varA, sep = " ")
    
    if(varA %in% colNames){
      augCoefMatrix[eq, varA] = coefA
    }
    
    sub = 1
    varB = paste("b", as.character(sub), sep = "")
    coefB = data[[1]][1]
    termB = paste(as.character(coefB), "*", varB, sep = " ")
    
    if(varB %in% colNames){
      augCoefMatrix[eq, varB] = coefB
    }
    
    sub = 1
    varC = paste("c", as.character(sub), sep = "")
    coefC = 1
    termC = paste(as.character(coefC), "*", varC, sep = " ")
    
    if(varC %in% colNames){
      augCoefMatrix[eq, varC] = coefC
    }
    
    sol = data[[2]][1]
    RHS[eq, 1] = sol
    
    equation = paste(termA, "+", termB, "+", termC, "=", sol)
    system[[equationCounter]] = equation
    equationCounter = equationCounter + 1
    
    eq = eq + 1
    
    # equation format 2
    sub = intervals
    varA = paste("a", as.character(sub), sep = "")
    coefA = data[[1]][dataPoints]^2
    termA = paste(as.character(coefA), "*", varA, sep = " ")
    
    if(varA %in% colNames){
      augCoefMatrix[eq, varA] = coefA
    }
    
    sub = intervals
    varB = paste("b", as.character(sub), sep = "")
    coefB = data[[1]][dataPoints]
    termB = paste(as.character(coefB), "*", varB, sep = " ")
    
    if(varB %in% colNames){
      augCoefMatrix[eq, varB] = coefB
    }
    
    sub = intervals
    varC = paste("c", as.character(sub), sep = "")
    coefC = 1
    termC = paste(as.character(coefC), "*", varC, sep = " ")
    
    if(varC %in% colNames){
      augCoefMatrix[eq, varC] = coefC
    }
    
    sol = data[[2]][dataPoints]
    RHS[eq, 1] = sol
    
    equation = paste(termA, "+", termB, "+", termC, "=", sol)
    system[[equationCounter]] = equation
    equationCounter = equationCounter + 1
    
    eq = eq + 1
    
    # condition 3
    i = 3 # 1 indexing
    while(i <= (dataPoints)){
      # equation format 1
      sub = i-1-1
      varA = paste("a", as.character(sub), sep = "")
      coefA = data[[1]][i-1]*2
      termA = paste(as.character(coefA), "*", varA, sep = " ")
      
      if(varA %in% colNames){
        augCoefMatrix[eq, varA] = coefA
      }
      
      sub = i-1-1
      varB = paste("b", as.character(sub), sep = "")
      coefB = 1
      termB = paste(as.character(coefB), "*", varB, sep = " ")
      
      if(varB %in% colNames){
        augCoefMatrix[eq, varB] = coefB
      }
      
      leftSide = paste(termA, "+", termB)
      
      sub = i-1
      varA = paste("a", as.character(sub), sep = "")
      coefA = data[[1]][i-1]*2
      termA = paste(as.character(coefA), "*", varA, sep = " ")
      
      if(varA %in% colNames){
        augCoefMatrix[eq, varA] = -(coefA)
      }
      
      sub = i-1
      varB = paste("b", as.character(sub), sep = "")
      coefB = 1
      termB = paste(as.character(coefB), "*", varB, sep = " ")
      
      if(varB %in% colNames){
        augCoefMatrix[eq, varB] = -(coefB)
      }
      
      rightSide = paste(termA, "+", termB)
      
      equation = paste(leftSide, "=", rightSide)
      system[[equationCounter]] = equation
      equationCounter = equationCounter + 1
      
      eq = eq + 1
      
      i = i + 1
    }
    
    sub = 1
    varA = paste("a", as.character(sub), sep = "")
    coefA = 1
    termA = paste(as.character(coefA), "*", varA, sep = " ")
    sol = 0
    
    lastEquation = paste(termA, "=", as.character(sol))
    system[[equationCounter]] = lastEquation
    
    finalMatrix <- cbind(augCoefMatrix, RHS)
    
    gaussJordanInitial <- matrix(c(0))
    gaussJordan <- matrix (gaussJordan(finalMatrix))
    gaussJordan = rbind(gaussJordanInitial, gaussJordan)
    
    i = 1
    j = i + 1
    while(!((data[[1]][i] <= x) & (data[[1]][j] >= x)) & (i <= intervals)){
      i = i + 1
      j = i + 1
    }
    
    gaussianStartingIndex = ((i - 1)*3) + 1
    z = 0
    exponent = 2
    summation = 0
    while(z < 3){
      p = gaussJordan[(gaussianStartingIndex + z),1]
      summation = summation + (p *(x^exponent))
      z = z + 1
      exponent = exponent - 1
    }
    
    polyPerInterval <- list(0)
    polyPerIntervalCounter = 1
    gaussCounter = 1
    while(gaussCounter <= length(gaussJordan)){
      qwerty = 1
      polynomial = ""
      while(qwerty <= 3){
        if(qwerty == 1){
          var = "x ^ 2"
        }else if(qwerty == 2){
          var = "x"
        }else{
          var = ""
        }
        
        if(qwerty == 3){
          term = gaussJordan[gaussCounter,1]
        }else{
          term = paste(gaussJordan[gaussCounter,1], "*", var, "+")
        }
        
        polynomial = paste(polynomial,term)
        qwerty = qwerty + 1
        gaussCounter = gaussCounter + 1
      }
      polyPerInterval[[polyPerIntervalCounter]] = polynomial
      polyPerIntervalCounter = polyPerIntervalCounter + 1
    }
    
    final <- list(
      # qsi.fxns = system,
      qsi.fxns = polyPerInterval,
      y = summation
    )
    
    return(final)
}

u <- shinyUI(pageWithSidebar(
  
  headerPanel("Quadratic Spline Interpolation"),
  sidebarPanel(
    p(em("First Vector")),
    textInput('vec1', 'Enter a vector (comma delimited)', "3.0, 4.5, 7.0, 9.0"),
    p(em("Second Vector")),
    textInput('vec2', 'Enter a vector (comma delimited)', "2.5, 1.0, 2.5, 0.5"),
    textInput('xValue', 'Enter a value for x', "5"),
    actionButton("go", "Calculate")
  ),
  
  
  
  mainPanel(
    h4('Answer'),
    verbatimTextOutput("qsi")
  )
))

s <- shinyServer(function(input, output) {

  calculate <- eventReactive(input$go, {
    vec1 <- c(as.numeric(unlist(strsplit(input$vec1,","))))
    vec2 <- c(as.numeric(unlist(strsplit(input$vec2,","))))
    data <- list(vec1,vec2)
    x <- as.numeric(input$xValue)
    if(length(vec1) > 1 & length(vec2) > 1){
      if(length(vec1) == length(vec2)){
        if(!is.unsorted(vec1)){
          if(x >= min(vec1) & x <= max(vec1)){
            poly.qsi <- poly.qsi(data, x)
            cat("=====> Polynomials\n\n")
            for(equations in poly.qsi[[1]]){
              cat(equations)
              cat("\n")
            }
            cat("\n=====> Estimated Value\n\n")
            cat(poly.qsi[[2]])
          }else{
            cat("ERROR!\nx Value must be inside the range of values of the First vector")
          }
        }else{
          cat("ERROR!\nFirst Vector must be in ASCENDING order")
        }
      }else{
        cat("ERROR!\nThe number of elements in the First and Second Vector must be EQUAL")
      }
    }else{
      cat("ERROR!\nThe number of elements in both Vectors must be GREATER THAN 1")
    }
    
  })
  
  output$qsi<-renderPrint({
    calculate()
  }
  )
  
}
)
shinyApp(ui = u, server = s)